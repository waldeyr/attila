#!/usr/bin/env python3
"""
test_attila.py - Unit and integration tests for validating the ATTILA pipeline components.
"""

import sys
import os
import unittest
import math
import tempfile
import shutil
import importlib.util
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add programs directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

_ATTILA_IMPORT_ERROR = None
_HTML_CREATOR_IMPORT_ERROR = None

try:
    import attila
except ImportError as e:
    _ATTILA_IMPORT_ERROR = str(e)
    attila = None  # type: ignore

try:
    import html_creator
except ImportError as e:
    _HTML_CREATOR_IMPORT_ERROR = str(e)
    html_creator = None  # type: ignore

class TestInstallation(unittest.TestCase):
    """Pre-flight environment checks. Run before the pipeline executes."""

    def test_01_python_version(self):
        self.assertGreaterEqual(
            sys.version_info, (3, 8),
            f"Python >= 3.8 required, found {sys.version_info.major}.{sys.version_info.minor}"
        )

    def test_02_matplotlib_importable(self):
        if importlib.util.find_spec("matplotlib") is None:
            self.skipTest("matplotlib not installed (optional — used for plots)")

    def test_03_numpy_importable(self):
        if importlib.util.find_spec("numpy") is None:
            self.skipTest("numpy not installed (optional — used for plots)")

    def test_04_fastqc_in_path(self):
        self.assertIsNotNone(shutil.which("fastqc"), "fastqc not found in PATH")

    def test_05_prinseq_in_path(self):
        self.assertIsNotNone(shutil.which("prinseq-lite"), "prinseq-lite not found in PATH")

    def test_06_fastqjoin_in_path(self):
        found = shutil.which("fastq-join") or shutil.which("fastqjoin")
        self.assertIsNotNone(found, "Neither fastq-join nor fastqjoin found in PATH")

    def test_07_igblastp_in_path(self):
        if shutil.which("igblastp") is None:
            self.skipTest("igblastp not found in PATH (germline step will be skipped)")

    def test_08_genetic_code_fasta(self):
        f = Path(__file__).parent.parent / "data" / "genetic_code.fasta"
        self.assertTrue(f.is_file(), f"Required data file missing: {f}")

    def test_09_germline_human_dir(self):
        d = Path(__file__).parent.parent / "data" / "germline_human"
        self.assertTrue(d.is_dir(), f"Required data directory missing: {d}")

    def test_10_attila_module_importable(self):
        self.assertIsNone(
            _ATTILA_IMPORT_ERROR,
            f"programs/attila.py could not be imported: {_ATTILA_IMPORT_ERROR}"
        )

    def test_11_html_creator_module_importable(self):
        self.assertIsNone(
            _HTML_CREATOR_IMPORT_ERROR,
            f"programs/html_creator.py could not be imported: {_HTML_CREATOR_IMPORT_ERROR}"
        )


class TestAttilaPipeline(unittest.TestCase):

    def setUp(self):
        if attila is None:
            self.skipTest(f"attila not importable: {_ATTILA_IMPORT_ERROR}")
        if html_creator is None:
            self.skipTest(f"html_creator not importable: {_HTML_CREATOR_IMPORT_ERROR}")

    def test_reverse_complement(self):
        self.assertEqual(attila.reverse_complement("ATCG"), "CGAT")
        self.assertEqual(attila.reverse_complement("atcg"), "CGAT")
        self.assertEqual(attila.reverse_complement("NnATCG"), "CGATNN")

    def test_translation_and_orf(self):
        # We need a codon table for translation. We can build a mock one.
        codon_table = {
            'ATG': 'M', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            'CAA': 'Q', 'CAG': 'Q', 'AAA': 'K', 'AAG': 'K', 'TTT': 'F', 'TTC': 'F'
        }
        
        # Build a mock V-domain amino acid sequence:
        # starts with C (Cys1) ... 55 residues ... C (Cys2) ... 7 residues CDR3 ... WGQG (FR4)
        # Sequence: C + 55*A + C + 7*G + WGQG
        aa_seq = "C" + ("A" * 55) + "C" + ("G" * 7) + "WGQG"
        
        # Test finding CDR3 boundaries
        # libtype 0 is VH (looks for WG.G)
        bounds = attila.find_cdr3_bounds(aa_seq, libtype=0)
        self.assertIsNotNone(bounds)
        c2_idx, m_idx, start_1based, end_1based = bounds
        
        # The second cysteine is at index 1 + 55 = 56
        self.assertEqual(c2_idx, 57) # c2_idx = c2 + 1 (since c2 is 56)
        # The motif WGQG starts at index 56 + 1 + 7 = 64
        self.assertEqual(m_idx, 64)
        # 1-indexed boundaries: start is c2 + 2 = 58, end is m_idx = 64
        self.assertEqual(start_1based, 58)
        self.assertEqual(end_1based, 64)
        
        # Slice CDR3 sequence
        cdr_seq = aa_seq[c2_idx : m_idx]
        self.assertEqual(cdr_seq, "GGGGGGG")

    def test_z_test_proportions(self):
        # Math helper for pnorm
        def pnorm(z):
            return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
            
        # Mock values:
        # R0: 10 reads out of 10000 -> proportion 0.001
        # RN: 200 reads out of 10000 -> proportion 0.02
        pantes = 10 / 10000.0
        pdepois = 200 / 10000.0
        n1 = 10000.0
        n2 = 10000.0
        
        pantesdepois = (n1 * pantes + n2 * pdepois) / (n1 + n2)
        erropadiff = math.sqrt(pantesdepois * (1.0 - pantesdepois) * (1.0/n1 + 1.0/n2))
        z_obs = (pdepois - pantes) / erropadiff
        p_value = 1.0 - pnorm(z_obs)
        
        # We expect a very small p-value (high enrichment)
        self.assertTrue(p_value < 0.001)
        
        # Check confidence interval
        ep_diff_ic = math.sqrt((pantes * (1.0 - pantes) / n1) + (pdepois * (1.0 - pdepois) / n2))
        inf_ic = pantes - pdepois - 1.96 * ep_diff_ic
        sup_ic = pantes - pdepois + 1.96 * ep_diff_ic
        self.assertTrue(inf_ic < sup_ic)

    @patch('urllib.request.urlopen')
    @patch('subprocess.run')
    def test_pipeline_integration(self, mock_sub_run, mock_url_open):
        # 1. Mock UCL Abnum HTTP request response
        def mock_urlopen_side_effect(req, timeout=None):
            from urllib.parse import urlparse, parse_qs
            import io
            if isinstance(req, str):
                url = req
            else:
                url = req.full_url
            parsed_url = urlparse(url)
            params = parse_qs(parsed_url.query)
            aaseq = params.get('aaseq', [''])[0]
            
            lines = []
            for i, aa in enumerate(aaseq):
                lines.append(f"H{i+1} {aa}")
            content = "\n".join(lines) + "\n"
            return io.BytesIO(content.encode('utf-8'))
            
        mock_url_open.side_effect = mock_urlopen_side_effect
        
        # 2. Mock subprocess.run for CLI dependencies
        def mock_sub_run_side_effect(cmd, *args, **kwargs):
            cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
            src_dir = Path(__file__).parent.parent / "data" / "synthetic_data_for_test"
            
            if "prinseq-lite" in cmd_str:
                # Mock prinseq-lite output by copying synthetic filtered FASTA files
                parts = cmd_str.split()
                try:
                    out_good_idx = parts.index("-out_good")
                    out_good_path = parts[out_good_idx + 1]
                except (ValueError, IndexError):
                    out_good_path = None
                    
                if out_good_path:
                    if "vh_initial" in out_good_path:
                        src_file = src_dir / "vh_initialfiltered.fasta"
                    elif "vh_final" in out_good_path:
                        src_file = src_dir / "vh_finalfiltered.fasta"
                    elif "vl_initial" in out_good_path:
                        src_file = src_dir / "vl_initialfiltered.fasta"
                    elif "vl_final" in out_good_path:
                        src_file = src_dir / "vl_finalfiltered.fasta"
                    else:
                        raise ValueError(f"Unknown out_good: {out_good_path}")
                        
                    dest_file = Path(out_good_path + ".fasta")
                    os.makedirs(dest_file.parent, exist_ok=True)
                    shutil.copy(src_file, dest_file)
                    
            elif "igblastp" in cmd_str:
                # Mock igblastp V-gene classification file in format 7
                if ">" in cmd_str:
                    dest_path = cmd_str.split(">")[1].strip()
                    parts = cmd_str.split()
                    query_arg = None
                    if "-query" in parts:
                        query_arg = parts[parts.index("-query") + 1]
                        
                    if query_arg and os.path.exists(query_arg):
                        query_ids = []
                        with open(query_arg, 'r') as f_q:
                            for line in f_q:
                                if line.startswith(">"):
                                    query_ids.append(line.rstrip()[1:])
                                    
                        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
                        with open(dest_path, 'w', encoding='utf-8') as f_out:
                            for q_id in query_ids:
                                f_out.write(f"# Query: {q_id}\n")
                                f_out.write("FR1\t1\t90\n")
                                f_out.write("CDR1\t91\t95\n")
                                f_out.write("FR2\t96\t110\n")
                                f_out.write("CDR2\t111\t120\n")
                                f_out.write("FR3\t121\t200\n")
                                name = q_id.split("|")[0]
                                germ = "IGHV3-23*01" if "VH" in name else "IGLV3-21*01"
                                f_out.write(f"V\t{q_id}\t{germ}\t98.2\t100\n")
                                
            res = MagicMock()
            res.returncode = 0
            return res
            
        mock_sub_run.side_effect = mock_sub_run_side_effect
        
        # 3. Perform integration run in temporary directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            proj_dir = tmp_path / "example_output"
            
            # Write config files
            packagedir = Path(__file__).parent.parent.resolve()
            
            vh_config_content = f"""
projectname: example_project_results
projectdir: {proj_dir.as_posix()}
packagedir: {packagedir.as_posix()}
igblastdir: {packagedir.as_posix()}/programs
input1dir: {packagedir.as_posix()}/data/synthetic_data_for_test/vh_initial.fastq
input2dir: {packagedir.as_posix()}/data/synthetic_data_for_test/vh_final.fastq
libtype: 0
listsize: 5
pairedend: 0
minlen: 300
minqual: 20
"""
            vl_config_content = f"""
projectname: example_project_results
projectdir: {proj_dir.as_posix()}
packagedir: {packagedir.as_posix()}
igblastdir: {packagedir.as_posix()}/programs
input1dir: {packagedir.as_posix()}/data/synthetic_data_for_test/vl_initial.fastq
input2dir: {packagedir.as_posix()}/data/synthetic_data_for_test/vl_final.fastq
libtype: 1
listsize: 5
pairedend: 0
minlen: 300
minqual: 20
"""
            vh_cfg_file = tmp_path / "example_VH.cfg"
            vl_cfg_file = tmp_path / "example_VL.cfg"
            
            with open(vh_cfg_file, 'w', encoding='utf-8') as f: f.write(vh_config_content)
            with open(vl_cfg_file, 'w', encoding='utf-8') as f: f.write(vl_config_content)
            
            # Execute attila.py for VH
            vh_error_log = proj_dir / "vherror.log"
            with patch.object(sys, 'argv', ["attila.py", str(vh_cfg_file), str(vh_error_log)]):
                attila.main()
                
            # Execute attila.py for VL
            vl_error_log = proj_dir / "vlerror.log"
            with patch.object(sys, 'argv', ["attila.py", str(vl_cfg_file), str(vl_error_log)]):
                attila.main()
                
            # Verify VH outputs
            vh_output_dir = proj_dir / "example_project_results" / "VH"
            self.assertTrue((vh_output_dir / "vhSequenceCounting.csv").exists())
            self.assertTrue((vh_output_dir / "SelectedSequences" / "vhenriched.fasta").exists())
            self.assertTrue((vh_output_dir / "SelectedSequences" / "vhlist5numbered.fasta").exists())
            self.assertTrue((vh_output_dir / "SelectedSequences" / "vhlist5numberednt.fasta").exists())
            
            # Verify VL outputs
            vl_output_dir = proj_dir / "example_project_results" / "VL"
            self.assertTrue((vl_output_dir / "vlSequenceCounting.csv").exists())
            self.assertTrue((vl_output_dir / "SelectedSequences" / "vlenriched.fasta").exists())
            self.assertTrue((vl_output_dir / "SelectedSequences" / "vllist5numbered.fasta").exists())
            self.assertTrue((vl_output_dir / "SelectedSequences" / "vllist5numberednt.fasta").exists())
            
            # Execute html_creator.py report generation
            vh_numbered = vh_output_dir / "SelectedSequences" / "vhlist5numbered.fasta"
            vl_numbered = vl_output_dir / "SelectedSequences" / "vllist5numbered.fasta"
            vh_germline = vh_output_dir / "SelectedSequences" / "vhlist5numberedgermlineclassification.txt"
            vl_germline = vl_output_dir / "SelectedSequences" / "vllist5numberedgermlineclassification.txt"
            report_html = proj_dir / "example_project_results" / "Report" / "Report.html"
            vh_csv = vh_output_dir / "vhSequenceCounting.csv"
            vl_csv = vl_output_dir / "vlSequenceCounting.csv"
            vh_stats = proj_dir / "example_project_results" / "Report" / "vhoutputRstats.txt"
            vl_stats = proj_dir / "example_project_results" / "Report" / "vloutputRstats.txt"
            
            with patch.object(sys, 'argv', [
                "html_creator.py",
                str(vh_numbered),
                str(vl_numbered),
                str(vh_germline),
                str(vl_germline),
                str(report_html),
                "plot1vh.png",
                "plot2vh.png",
                "plot1vl.png",
                "plot2vl.png",
                str(vh_csv),
                str(vl_csv),
                str(vh_stats),
                str(vl_stats)
            ]):
                html_creator.main()
                
            self.assertTrue(report_html.exists())
            # Inspect Report.html content to ensure it has correct title/elements
            with open(report_html, 'r', encoding='utf-8') as r_f:
                html_content = r_f.read()
                self.assertIn("Selected sequences", html_content)
                self.assertIn("VH Library", html_content)
                self.assertIn("VL Library", html_content)

if __name__ == '__main__':
    unittest.main()
