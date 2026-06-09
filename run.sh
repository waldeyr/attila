#!/bin/bash

# =====================================================================
# ATTILA - AutomaTed Tool For Immunoglobulin Analysis
# Wrapper & Orchestrator Script
# =====================================================================

set -e

# Core paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROGRAMS_DIR="${SCRIPT_DIR}/programs"
DATA_DIR="${SCRIPT_DIR}/data"

# Add igblastp to PATH if installed but not yet exported (e.g. fresh terminal after install.sh)
if ! command -v igblastp &>/dev/null; then
    for _igblast_dir in /opt/ncbi-igblast-*/bin; do
        if [ -x "${_igblast_dir}/igblastp" ]; then
            export PATH="${_igblast_dir}:${PATH}"
            break
        fi
    done
fi

# Default settings
CONFIG_FILE=""
STEPS="all"
INTERACTIVE=false
LIB_TYPE="" # vh, vl, both, or empty (read from config)
RUN_EXAMPLE=false

show_help() {
    echo "====================================================================="
    echo "  ATTILA: Automated Tool for Immunoglobulin Analysis"
    echo "====================================================================="
    echo "Usage: ./run.sh [options]"
    echo ""
    echo "Options:"
    echo "  -c, --config FILE       Path to configuration file (.cfg)"
    echo "  -s, --steps STEPS       Comma-separated steps to execute"
    echo "                          Available steps: join, filter, translate, frequency, number, nt-recovery, stats, germline, report"
    echo "                          Default: all"
    echo "  -t, --type TYPE         Library type to run: vh, vl, or both"
    echo "  -i, --interactive       Run in interactive configuration wizard mode"
    echo "  -e, --example           Run the pipeline using synthetic demo data"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Verify installation using built-in synthetic data (no config needed)"
    echo "  ./run.sh --example"
    echo "  # Output is saved to: example_output/example_project_results/"
    echo ""
    echo "  # Run full pipeline with config file"
    echo "  ./run.sh --config my_project_VH.cfg"
    echo ""
    echo "  # Run both VH and VL configs sequentially"
    echo "  ./run.sh --config my_project_VH.cfg --type both"
    echo ""
    echo "  # Run only filter and translation steps"
    echo "  ./run.sh --config my_project_VH.cfg --steps filter,translate"
    echo ""
    echo "  # Run interactive wizard"
    echo "  ./run.sh --interactive"
    echo "====================================================================="
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -c|--config)
            CONFIG_FILE="$2"
            shift; shift
            ;;
        -s|--steps)
            STEPS="$2"
            shift; shift
            ;;
        --all)
            STEPS="all"
            shift
            ;;
        -t|--type)
            LIB_TYPE=$(echo "$2" | tr '[:upper:]' '[:lower:]')
            shift; shift
            ;;
        -i|--interactive)
            INTERACTIVE=true
            shift
            ;;
        -e|--example)
            RUN_EXAMPLE=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            show_help
            exit 1
            ;;
    esac
done

parse_config_value() {
    local key="$1"
    local file="$2"
    grep "^${key}:" "$file" | sed "s/^${key}:[[:space:]]*//g" | sed 's/[[:space:]]*$//'
}

check_dependencies() {
    echo "--------------------------------------------------------"
    echo "Checking CLI dependencies..."
    echo "--------------------------------------------------------"
    local missing=()
    
    for tool in fastqc prinseq-lite python3; do
        if ! command -v "$tool" &> /dev/null; then
            missing+=("$tool")
        fi
    done
    
    if ! command -v fastq-join &> /dev/null && ! command -v fastqjoin &> /dev/null; then
        missing+=("fastq-join")
    fi
    
    if [ ${#missing[@]} -ne 0 ]; then
        echo "WARNING: The following dependencies are missing from your PATH:"
        for m in "${missing[@]}"; do
            echo "  - $m"
        done
        echo "Please make sure they are installed for the pipeline to complete."
        echo ""
    else
        echo "All core CLI dependencies found!"
    fi
}

run_preflight_checks() {
    echo "--------------------------------------------------------"
    echo "Running installation pre-flight checks..."
    echo "--------------------------------------------------------"
    set +e
    python3 "${PROGRAMS_DIR}/test_attila.py" TestInstallation -v 2>&1
    local preflight_exit=$?
    set -e
    if [ $preflight_exit -ne 0 ]; then
        echo ""
        echo "ERROR: Pre-flight installation checks FAILED."
        echo "Run the installer to fix missing dependencies:"
        echo "  ./install.sh"
        echo ""
        echo "To run the full test suite manually: python3 programs/test_attila.py -v"
        exit $preflight_exit
    fi
    echo "Pre-flight checks passed."
    echo "--------------------------------------------------------"
}

run_interactive_wizard() {
    echo "********************************************************************************"
    echo "             ATTILA Interactive Configuration Wizard"
    echo "********************************************************************************"
    
    PACKAGEDIR="${SCRIPT_DIR}/"
    
    read -p "Enter project name: " PROJECT_NAME
    
    while true; do
        read -p "Enter project path (where results will be saved): " PROJECT_DIR
        PROJECT_DIR="${PROJECT_DIR/#\~/$HOME}"
        if [[ -d "$PROJECT_DIR" ]]; then
            break
        else
            read -p "Directory does not exist. Create it? (y/n, default: y): " mkdir_confirm
            mkdir_confirm=${mkdir_confirm:-y}
            if [[ "$mkdir_confirm" =~ ^[Yy]$ ]]; then
                mkdir -p "$PROJECT_DIR"
                break
            fi
        fi
    done
    
    read -p "Reads are paired-end? (y/n, default: y): " PAIRED_YN
    PAIRED_YN=${PAIRED_YN:-y}
    if [[ "$PAIRED_YN" =~ ^[Nn]$ ]]; then
        PAIRED=0
    else
        PAIRED=1
    fi
    
    declare -A FILE_PATHS
    
    if [ $PAIRED -eq 1 ]; then
        for v in VH VL; do
            for r in R0 RN; do
                for read_idx in 1 2; do
                    while true; do
                        read -p "Enter path for $v $r read r$read_idx (.fastq/.fq): " file_path
                        file_path="${file_path/#\~/$HOME}"
                        if [[ -f "$file_path" ]]; then
                            FILE_PATHS["${v}_${r}_r${read_idx}"]="$file_path"
                            break
                        else
                            echo "Error: File does not exist. Please enter a valid path."
                        fi
                    done
                done
            done
        done
    else
        for v in VH VL; do
            for r in R0 RN; do
                while true; do
                    read -p "Enter path for $v $r (.fastq/.fq): " file_path
                    file_path="${file_path/#\~/$HOME}"
                    if [[ -f "$file_path" ]]; then
                        FILE_PATHS["${v}_${r}"]="$file_path"
                        break
                    else
                        echo "Error: File does not exist. Please enter a valid path."
                    fi
                done
            done
        done
    fi
    
    while true; do
        read -p "Enter IgBlast package path: " IGBLAST_DIR
        IGBLAST_DIR="${IGBLAST_DIR/#\~/$HOME}"
        if [[ -d "$IGBLAST_DIR" ]]; then
            break
        else
            echo "Error: Directory does not exist. Please enter a valid path."
        fi
    done
    
    read -p "Minimum read length (default: 300): " MIN_LEN
    MIN_LEN=${MIN_LEN:-300}
    
    read -p "Minimum base quality (default: 20): " MIN_QUAL
    MIN_QUAL=${MIN_QUAL:-20}
    
    while true; do
        read -p "Enter number of candidates to rank (integer): " CANDIDATES
        if [[ "$CANDIDATES" =~ ^[0-9]+$ ]]; then
            break
        else
            echo "Error: Please enter a positive integer."
        fi
    done
    
    VH_CFG="${PROJECT_DIR}/${PROJECT_NAME}_VH.cfg"
    VL_CFG="${PROJECT_DIR}/${PROJECT_NAME}_VL.cfg"
    
    echo "Writing configuration files..."
    for v in VH VL; do
        CFG_FILE=""
        if [ "$v" == "VH" ]; then
            CFG_FILE="$VH_CFG"
            libtype=0
        else
            CFG_FILE="$VL_CFG"
            libtype=1
        fi
        
        {
            echo "----------------------------------------------------------------------"
            echo "# Settings for immunoglobulin sequence analysis"
            echo "# [ Section: files and directories ]"
            echo "projectname: ${PROJECT_NAME}"
            echo "projectdir: ${PROJECT_DIR}"
            echo "packagedir: ${PACKAGEDIR}"
            echo "igblastdir: ${IGBLAST_DIR}"
            
            if [ $PAIRED -eq 1 ]; then
                echo "input1r1dir: ${FILE_PATHS[${v}_R0_r1]}"
                echo "input1r2dir: ${FILE_PATHS[${v}_R0_r2]}"
                echo "input2r1dir: ${FILE_PATHS[${v}_RN_r1]}"
                echo "input2r2dir: ${FILE_PATHS[${v}_RN_r2]}"
            else
                echo "input1dir: ${FILE_PATHS[${v}_R0]}"
                echo "input2dir: ${FILE_PATHS[${v}_RN]}"
            fi
            
            echo "# [ Section: analysis arguments ]"
            echo "libtype: ${libtype}"
            echo "listsize: ${CANDIDATES}"
            echo "pairedend: ${PAIRED}"
            echo "minlen: ${MIN_LEN}"
            echo "minqual: ${MIN_QUAL}"
        } > "$CFG_FILE"
        echo "Created: $CFG_FILE"
    done
    
    read -p "Configuration complete. Run the analysis now? (y/n, default: y): " run_now
    run_now=${run_now:-y}
    if [[ ! "$run_now" =~ ^[Nn]$ ]]; then
        CONFIG_FILE="$VH_CFG"
        LIB_TYPE="both"
    else
        echo "Configuration files written to ${PROJECT_DIR}. You can run them later using:"
        echo "  ./run.sh --config ${VH_CFG}"
        echo "  ./run.sh --config ${VL_CFG}"
        exit 0
    fi
}

compile_report() {
    local proj_dir="$1"
    local proj_name="$2"
    local report_dir="${proj_dir}/${proj_name}/Report"
    
    echo "--------------------------------------------------------"
    echo "Compiling final HTML report..."
    echo "--------------------------------------------------------"
    
    local vh_numbered=$(find "${proj_dir}/${proj_name}/VH" -name "*numbered.fasta" 2>/dev/null | head -n 1)
    local vh_germline=$(find "${proj_dir}/${proj_name}/VH" -name "*germlineclassification.txt" 2>/dev/null | head -n 1)
    local vl_numbered=$(find "${proj_dir}/${proj_name}/VL" -name "*numbered.fasta" 2>/dev/null | head -n 1)
    local vl_germline=$(find "${proj_dir}/${proj_name}/VL" -name "*germlineclassification.txt" 2>/dev/null | head -n 1)
    
    local plot1vh="plot1vh.r"
    local plot2vh="plot2vh.r"
    local plot1vl="plot1vl.r"
    local plot2vl="plot2vl.r"
    
    if [[ -f "$vh_numbered" && -f "$vl_numbered" && -f "$vh_germline" && -f "$vl_germline" ]]; then
        echo "Found VH and VL files. Generating report..."
        python3 "${PROGRAMS_DIR}/html_creator.py" \
            "$vh_numbered" \
            "$vl_numbered" \
            "$vh_germline" \
            "$vl_germline" \
            "${proj_dir}/${proj_name}/index.html" \
            "Report/plot1vh.png" \
            "Report/plot2vh.png" \
            "Report/plot1vl.png" \
            "Report/plot2vl.png" \
            "${proj_dir}/${proj_name}/VH/vhSequenceCounting.csv" \
            "${proj_dir}/${proj_name}/VL/vlSequenceCounting.csv" \
            "${report_dir}/vhoutputRstats.txt" \
            "${report_dir}/vloutputRstats.txt" > "${report_dir}/webpage.log" 2>&1
            
        if [ $? -eq 0 ]; then
            echo "Report generated successfully: ${report_dir}/Report.html"
        else
            echo "Error generating Report.html. Check log: ${report_dir}/webpage.log"
        fi
    else
        echo "Warning: Could not compile final HTML report because some VH or VL output files are missing."
        echo "Make sure both VH and VL analyses have completed successfully."
    fi
}

run_pipeline() {
    local cfg="$1"
    
    if [[ ! -f "$cfg" ]]; then
        echo "Error: Configuration file not found at '$cfg'"
        exit 1
    fi
    
    local proj_name=$(parse_config_value "projectname" "$cfg")
    local proj_dir=$(parse_config_value "projectdir" "$cfg")
    local pack_dir=$(parse_config_value "packagedir" "$cfg")
    local lib_val=$(parse_config_value "libtype" "$cfg")
    
    local lib_name="VH"
    if [ "$lib_val" -eq 1 ]; then
        lib_name="VL"
    fi
    
    echo "========================================================"
    echo "Running ATTILA pipeline for $lib_name..."
    echo "Project: $proj_name"
    echo "Destination: $proj_dir"
    echo "Steps to run: $STEPS"
    echo "========================================================"
    
    # 1. Ensure target directories exist
    local full_proj_dir="${proj_dir}/${proj_name}"
    mkdir -p "${full_proj_dir}/${lib_name}"
    mkdir -p "${full_proj_dir}/Report"
    
    # 2. Copy generic logos/images from img folder to Report
    cp "${SCRIPT_DIR}/img/attila.png" "${full_proj_dir}/Report/app5.png" || true
    cp "${SCRIPT_DIR}/img/bioinfo.png" "${full_proj_dir}/Report/bioinfo.png" || true
    cp "${SCRIPT_DIR}/img/unb.png" "${full_proj_dir}/Report/unb_logo1.png" || true
    
    # 3. Setup environment steps variable
    export ATTILA_STEPS="$STEPS"
    
    # 4. Execute Python pipeline runner
    local error_log="${full_proj_dir}/$(echo "$lib_name" | tr '[:upper:]' '[:lower:]')error.log"
    
    echo "Launching attila.py..."
    local start_time=$(date +%s)
    
    # Run Python script
    set +e
    python3 "${PROGRAMS_DIR}/attila.py" "$cfg" "$error_log"
    local exit_code=$?
    set -e
    
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    
    if [ $exit_code -eq 0 ]; then
        echo "--------------------------------------------------------"
        echo "$lib_name Analysis Completed Successfully in ${elapsed}s!"
        echo "--------------------------------------------------------"
    else
        echo "--------------------------------------------------------"
        echo "ERROR: $lib_name Analysis failed with exit code $exit_code."
        echo "Check log file: $error_log"
        echo "--------------------------------------------------------"
        exit $exit_code
    fi
}

run_example() {
    echo "========================================================"
    echo "Running ATTILA demonstration with synthetic data..."
    echo "========================================================"
    
    local workspace_dir="${SCRIPT_DIR}"
    local example_dir="${workspace_dir}/example_output"
    
    # Check if dependencies are present
    local has_qc=true
    for tool in fastqc prinseq-lite; do
        if ! command -v "$tool" &>/dev/null; then
            has_qc=false
        fi
    done
    
    rm -rf "${example_dir}"
    mkdir -p "${example_dir}"
    
    # Generate configurations from templates
    echo "Generating example configuration files..."
    sed "s|@WORKSPACE_DIR@|${workspace_dir}|g" "${workspace_dir}/data/synthetic_data_for_test/example_VH.cfg" > "${example_dir}/example_VH.cfg"
    sed "s|@WORKSPACE_DIR@|${workspace_dir}|g" "${workspace_dir}/data/synthetic_data_for_test/example_VL.cfg" > "${example_dir}/example_VL.cfg"
    
    if [ "$has_qc" = true ]; then
        echo "Dependencies found. Running FULL pipeline with synthetic FASTQ reads..."
        local steps_to_run="join,filter,translate,frequency,number,nt-recovery,stats,report"
        if command -v igblastp &>/dev/null; then
            steps_to_run="${steps_to_run},germline"
        fi
        export STEPS="${steps_to_run}"
    else
        echo "WARNING: fastqc and/or prinseq-lite not found in PATH."
        echo "Running in FALLBACK mode: executing Python-only stages using pre-filtered FASTA data..."
        
        # Setup output folders
        mkdir -p "${example_dir}/example_project_results/VH/InitialRound"
        mkdir -p "${example_dir}/example_project_results/VH/FinalRound"
        mkdir -p "${example_dir}/example_project_results/VL/InitialRound"
        mkdir -p "${example_dir}/example_project_results/VL/FinalRound"
        mkdir -p "${example_dir}/example_project_results/Report"
        
        # Copy pre-filtered FASTA files
        cp "${workspace_dir}/data/synthetic_data_for_test/vh_initialfiltered.fasta" "${example_dir}/example_project_results/VH/InitialRound/vh_initialfiltered.fasta"
        cp "${workspace_dir}/data/synthetic_data_for_test/vh_finalfiltered.fasta" "${example_dir}/example_project_results/VH/FinalRound/vh_finalfiltered.fasta"
        cp "${workspace_dir}/data/synthetic_data_for_test/vl_initialfiltered.fasta" "${example_dir}/example_project_results/VL/InitialRound/vl_initialfiltered.fasta"
        cp "${workspace_dir}/data/synthetic_data_for_test/vl_finalfiltered.fasta" "${example_dir}/example_project_results/VL/FinalRound/vl_finalfiltered.fasta"
        
        export STEPS="translate,frequency,number,nt-recovery,stats,report"
    fi
    
    # Touch empty germline files if germline step is not run, to allow compiling report
    if [[ ! "$STEPS" == *germline* ]]; then
        mkdir -p "${example_dir}/example_project_results/VH/SelectedSequences"
        mkdir -p "${example_dir}/example_project_results/VL/SelectedSequences"
        touch "${example_dir}/example_project_results/VH/SelectedSequences/vhlist5numberedgermlineclassification.txt"
        touch "${example_dir}/example_project_results/VL/SelectedSequences/vllist5numberedgermlineclassification.txt"
    fi
    
    CONFIG_FILE="${example_dir}/example_VH.cfg"
    LIB_TYPE="both"
}

main() {
    # Pre-flight: verify Python environment, data files, and CLI tools
    run_preflight_checks

    # Check dependencies (warns only, does not fail)
    check_dependencies
    
    if [ "$RUN_EXAMPLE" = true ]; then
        run_example
    elif [ "$INTERACTIVE" = true ]; then
        run_interactive_wizard
    fi
    
    if [[ -z "$CONFIG_FILE" ]]; then
        # If no config file is passed, fall back to interactive wizard
        echo "No configuration file specified."
        run_interactive_wizard
    fi
    
    # Resolve absolute path for config file
    CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")
    
    # Determine which config files to run
    local configs_to_run=()
    
    if [ "$LIB_TYPE" == "both" ]; then
        if [[ "$CONFIG_FILE" == *VH.cfg ]]; then
            configs_to_run+=("$CONFIG_FILE")
            configs_to_run+=("${CONFIG_FILE/VH.cfg/VL.cfg}")
        elif [[ "$CONFIG_FILE" == *VL.cfg ]]; then
            configs_to_run+=("${CONFIG_FILE/VL.cfg/VH.cfg}")
            configs_to_run+=("$CONFIG_FILE")
        else
            # Try to search or run it as is, and warn
            echo "Warning: Running single config since type was both but naming did not match VH/VL pattern."
            configs_to_run+=("$CONFIG_FILE")
        fi
    else
        configs_to_run+=("$CONFIG_FILE")
    fi
    
    # Execute the pipelines
    for cfg in "${configs_to_run[@]}"; do
        run_pipeline "$cfg"
    done
    
    # Post-execution Report Compilation
    if [[ "$STEPS" == "all" || "$STEPS" == *report* ]]; then
        local first_cfg="${configs_to_run[0]}"
        local proj_name=$(parse_config_value "projectname" "$first_cfg")
        local proj_dir=$(parse_config_value "projectdir" "$first_cfg")
        compile_report "$proj_dir" "$proj_name"
    fi
}

main
