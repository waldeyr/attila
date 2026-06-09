#!/bin/bash

# =====================================================================
# ATTILA - Dependency Installer
# Detects macOS or Linux (Debian/Ubuntu via apt, RHEL/Fedora/Rocky via dnf)
# and installs all required tools:
#   fastqc, fastq-join (ea-utils), prinseq-lite, igblastp
#   Python packages: matplotlib, numpy
# =====================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IGBLAST_VERSION="1.22.0"
PRINSEQ_VERSION="0.20.4"

# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------

info()    { echo "[INFO]  $*"; }
success() { echo "[OK]    $*"; }
warn()    { echo "[WARN]  $*"; }
error()   { echo "[ERROR] $*" >&2; exit 1; }

already_installed() {
    local tool="$1"
    if command -v "$tool" &>/dev/null; then
        success "$tool already installed: $(command -v "$tool")"
        return 0
    fi
    return 1
}

# ----------------------------------------------------------------
# Detect OS and Linux distro family
# ----------------------------------------------------------------

detect_os() {
    case "$(uname -s)" in
        Darwin)  echo "macos" ;;
        Linux)   echo "linux" ;;
        *)       error "Unsupported OS: $(uname -s)" ;;
    esac
}

detect_arch() {
    case "$(uname -m)" in
        x86_64|amd64) echo "x64" ;;
        arm64|aarch64) echo "arm64" ;;
        *) echo "x64" ;;  # fallback
    esac
}

# Returns: apt | dnf | unknown
detect_linux_pkg_manager() {
    if command -v apt-get &>/dev/null; then
        echo "apt"
    elif command -v dnf &>/dev/null; then
        echo "dnf"
    elif command -v yum &>/dev/null; then
        echo "dnf"  # treat yum as dnf-compatible (same package names)
    else
        echo "unknown"
    fi
}

OS=$(detect_os)
ARCH=$(detect_arch)

if [ "$OS" = "linux" ]; then
    LINUX_PKG=$(detect_linux_pkg_manager)
    info "Detected OS: Linux — package manager: $LINUX_PKG ($ARCH)"
else
    LINUX_PKG=""
    info "Detected OS: $OS ($ARCH)"
fi
echo ""

# ----------------------------------------------------------------
# macOS: require Homebrew
# ----------------------------------------------------------------

check_brew() {
    if ! command -v brew &>/dev/null; then
        error "Homebrew is required on macOS. Install it from https://brew.sh and re-run this script."
    fi
}

# ----------------------------------------------------------------
# Linux: package manager helpers
# ----------------------------------------------------------------

pkg_install() {
    # Usage: pkg_install pkg1 pkg2 ...
    case "$LINUX_PKG" in
        apt) sudo apt-get install -y "$@" ;;
        dnf)
            if command -v dnf &>/dev/null; then
                sudo dnf install -y "$@"
            else
                sudo yum install -y "$@"
            fi
            ;;
        *) error "No supported package manager found (apt/dnf/yum)." ;;
    esac
}

# ----------------------------------------------------------------
# Python packages (both OS)
# ----------------------------------------------------------------

install_python_packages() {
    info "Installing Python packages (matplotlib, numpy)..."
    if python3 -m pip install --upgrade matplotlib numpy; then
        success "Python packages installed."
    else
        warn "pip install failed. Try manually: pip3 install matplotlib numpy"
    fi
}

# ----------------------------------------------------------------
# prinseq-lite (Perl script, both OS)
# ----------------------------------------------------------------

install_prinseq() {
    if already_installed prinseq-lite; then return; fi

    info "Installing prinseq-lite ${PRINSEQ_VERSION}..."

    local url="https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-${PRINSEQ_VERSION}.tar.gz"
    local tmp_dir
    tmp_dir=$(mktemp -d)

    if ! command -v perl &>/dev/null; then
        if [ "$OS" = "macos" ]; then
            info "Installing perl via Homebrew..."
            brew install perl
        else
            info "Installing perl..."
            pkg_install perl
        fi
    fi

    info "Downloading prinseq-lite from SourceForge..."
    wget --no-check-certificate "$url" -O "${tmp_dir}/prinseq.tar.gz" \
        || curl -L "$url" -o "${tmp_dir}/prinseq.tar.gz" \
        || error "Failed to download prinseq-lite. Check your internet connection."

    tar -xzf "${tmp_dir}/prinseq.tar.gz" -C "${tmp_dir}"
    sudo mv "${tmp_dir}/prinseq-lite-${PRINSEQ_VERSION}/prinseq-lite.pl" /usr/local/bin/prinseq-lite
    sudo chmod +x /usr/local/bin/prinseq-lite
    rm -rf "${tmp_dir}"

    success "prinseq-lite installed at /usr/local/bin/prinseq-lite"
}

# ----------------------------------------------------------------
# igblastp (optional — germline step)
# ----------------------------------------------------------------

install_igblastp() {
    if already_installed igblastp; then return; fi

    info "Installing NCBI IgBlast ${IGBLAST_VERSION} (optional — for germline classification)..."

    local base_url="ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}"
    local tmp_dir
    tmp_dir=$(mktemp -d)
    local install_dir="/opt/ncbi-igblast-${IGBLAST_VERSION}"
    local tarball

    if [ "$OS" = "macos" ]; then
        if [ "$ARCH" = "arm64" ]; then
            tarball="ncbi-igblast-${IGBLAST_VERSION}-x64-macosx.tar.gz"
            warn "No native ARM64 IgBlast binary available. Downloading x64 binary (requires Rosetta 2)."
            if ! /usr/bin/pgrep -q oahd 2>/dev/null; then
                warn "Rosetta 2 may not be installed. Run: softwareupdate --install-rosetta --agree-to-license"
            fi
        else
            tarball="ncbi-igblast-${IGBLAST_VERSION}-x64-macosx.tar.gz"
        fi
    else
        tarball="ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz"
        # Install libgomp if missing (required by igblastp on Linux)
        if ! ldconfig -p 2>/dev/null | grep -q libgomp; then
            info "Installing libgomp..."
            case "$LINUX_PKG" in
                apt) pkg_install libgomp1 ;;
                dnf) pkg_install libgomp ;;
            esac
        fi
    fi

    info "Downloading IgBlast: ${tarball}..."
    wget "${base_url}/${tarball}" -O "${tmp_dir}/igblast.tar.gz" \
        || curl -L "${base_url}/${tarball}" -o "${tmp_dir}/igblast.tar.gz" \
        || { warn "Failed to download IgBlast. Skipping (germline step will be unavailable)."; rm -rf "${tmp_dir}"; return; }

    sudo mkdir -p /opt
    sudo tar -xzf "${tmp_dir}/igblast.tar.gz" -C /opt
    sudo chmod +x "${install_dir}/bin/"*
    rm -rf "${tmp_dir}"

    # Add to PATH via a profile snippet
    local profile_snippet="export PATH=\"${install_dir}/bin:\$PATH\""
    for rc_file in "${HOME}/.bashrc" "${HOME}/.zshrc"; do
        if [ -f "$rc_file" ] && ! grep -q "ncbi-igblast" "$rc_file"; then
            echo "" >> "$rc_file"
            echo "# IgBlast (added by ATTILA installer)" >> "$rc_file"
            echo "$profile_snippet" >> "$rc_file"
            info "Added IgBlast to PATH in $rc_file"
        fi
    done

    export PATH="${install_dir}/bin:${PATH}"
    success "igblastp installed at ${install_dir}/bin/igblastp"
    warn "Restart your shell (or run: export PATH=\"${install_dir}/bin:\$PATH\") to use igblastp in new sessions."
}

# ----------------------------------------------------------------
# macOS installation
# ----------------------------------------------------------------

install_macos() {
    check_brew

    echo "========================================================"
    echo "Installing ATTILA dependencies on macOS via Homebrew"
    echo "========================================================"
    echo ""

    # fastqc
    if already_installed fastqc; then :
    else
        info "Installing fastqc..."
        brew install fastqc
        success "fastqc installed."
    fi

    # fastq-join (ea-utils) — not in Homebrew; compile from source
    if already_installed fastq-join || already_installed fastqjoin; then :
    else
        info "fastq-join is not available via Homebrew. Building from source (ea-utils)..."

        # Ensure build tools available
        if ! command -v gcc &>/dev/null; then
            info "Installing Xcode command-line tools..."
            xcode-select --install 2>/dev/null || true
            warn "If prompted, complete the Xcode tools installation and re-run this script."
        fi

        local tmp_dir
        tmp_dir=$(mktemp -d)

        info "Cloning ea-utils from GitHub..."
        git clone https://github.com/ExpressionAnalysis/ea-utils.git "${tmp_dir}/ea-utils" \
            || error "Failed to clone ea-utils. Check your internet connection and git installation."

        (
            cd "${tmp_dir}/ea-utils/clipper"
            make fastq-join 2>/dev/null || make 2>/dev/null
        )

        if [ -f "${tmp_dir}/ea-utils/clipper/fastq-join" ]; then
            sudo cp "${tmp_dir}/ea-utils/clipper/fastq-join" /usr/local/bin/fastq-join
            sudo chmod +x /usr/local/bin/fastq-join
            success "fastq-join installed at /usr/local/bin/fastq-join"
        else
            warn "Build failed. Trying pre-compiled binary approach..."
            # Fallback: homebrew tap if ever added, or warn user
            warn "Could not build fastq-join automatically."
            warn "Manual option: conda install -c bioconda ea-utils"
            warn "              OR: docker run ... (see README)"
        fi
        rm -rf "${tmp_dir}"
    fi

    # prinseq-lite
    install_prinseq

    # Python packages
    install_python_packages

    # igblastp (optional)
    install_igblastp
}

# ----------------------------------------------------------------
# Linux installation (apt: Debian/Ubuntu | dnf: RHEL/Fedora/Rocky)
# ----------------------------------------------------------------

install_linux() {
    if [ "$LINUX_PKG" = "unknown" ]; then
        error "No supported package manager found. Supported: apt (Debian/Ubuntu), dnf/yum (RHEL/Fedora/CentOS/Rocky/Alma)."
    fi

    echo "========================================================"
    echo "Installing ATTILA dependencies on Linux ($LINUX_PKG)"
    echo "========================================================"
    echo ""

    # Update package lists
    case "$LINUX_PKG" in
        apt)
            info "Updating package lists..."
            sudo apt-get update -qq
            ;;
        dnf)
            info "Enabling EPEL repository (required for fastqc on RHEL-based systems)..."
            # EPEL provides fastqc on RHEL/CentOS/Rocky/Alma
            if command -v dnf &>/dev/null; then
                sudo dnf install -y epel-release 2>/dev/null || true
            else
                sudo yum install -y epel-release 2>/dev/null || true
            fi
            ;;
    esac

    # fastqc
    # apt: 'fastqc'  |  dnf: 'fastqc' (via EPEL)
    if already_installed fastqc; then :
    else
        info "Installing fastqc..."
        pkg_install fastqc
        success "fastqc installed."
    fi

    # fastq-join
    # apt: 'ea-utils'  |  dnf: 'ea-utils' (via EPEL) — if not available, build from source
    if already_installed fastq-join || already_installed fastqjoin; then :
    else
        info "Installing ea-utils (provides fastq-join)..."
        if pkg_install ea-utils 2>/dev/null; then
            success "fastq-join installed."
        else
            warn "ea-utils not found in repositories. Building fastq-join from source..."
            _build_fastqjoin_from_source
        fi
    fi

    # wget (used by prinseq and igblast downloads)
    if ! command -v wget &>/dev/null; then
        pkg_install wget
    fi

    # prinseq-lite (Perl script — same method on all distros)
    install_prinseq

    # Python packages
    install_python_packages

    # igblastp (optional)
    install_igblastp
}

_build_fastqjoin_from_source() {
    local tmp_dir
    tmp_dir=$(mktemp -d)

    # Ensure build tools are present
    case "$LINUX_PKG" in
        apt) pkg_install gcc make git ;;
        dnf) pkg_install gcc make git ;;
    esac

    info "Cloning ea-utils from GitHub..."
    git clone https://github.com/ExpressionAnalysis/ea-utils.git "${tmp_dir}/ea-utils" \
        || error "Failed to clone ea-utils. Check your internet connection."

    (
        cd "${tmp_dir}/ea-utils/clipper"
        make fastq-join 2>/dev/null || make 2>/dev/null
    )

    if [ -f "${tmp_dir}/ea-utils/clipper/fastq-join" ]; then
        sudo cp "${tmp_dir}/ea-utils/clipper/fastq-join" /usr/local/bin/fastq-join
        sudo chmod +x /usr/local/bin/fastq-join
        success "fastq-join built and installed at /usr/local/bin/fastq-join"
    else
        warn "Build failed. fastq-join could not be installed automatically."
        warn "Manual option: conda install -c bioconda ea-utils"
        warn "              OR use Docker: see README"
    fi
    rm -rf "${tmp_dir}"
}

# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------

main() {
    echo "====================================================================="
    echo "  ATTILA Dependency Installer"
    echo "====================================================================="
    echo ""

    case "$OS" in
        macos) install_macos ;;
        linux) install_linux ;;
    esac

    echo ""
    echo "====================================================================="
    echo "  Installation complete. Verifying with pre-flight checks..."
    echo "====================================================================="
    echo ""
    python3 "${SCRIPT_DIR}/programs/test_attila.py" TestInstallation -v 2>&1 || true
    echo ""
    info "If any checks failed, review the messages above."
    info "To run the full ATTILA test suite: python3 programs/test_attila.py -v"
    info "To start the pipeline:             ./run.sh --example"
}

main
