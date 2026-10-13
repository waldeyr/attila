"""Known defects, documented as strict xfails.

Each test asserts the *correct* behaviour, so it fails today. When a bug is
fixed the test starts passing, xfail_strict turns that into a failure, and the
marker should be removed.
"""

import pytest

import attila
from helpers import REPO_ROOT, lib_dirs, write


@pytest.mark.xfail(raises=FileNotFoundError, reason=(
    "paired-end: step_filter names its output after Path('X.fqjoin').stem = 'X' (attila.py L242), "
    "but step_translate reads 'X.fqjoinfiltered.fasta' (L265-271)"))
def test_paired_end_filter_output_feeds_translate(make_cfg, fake_tools, tmp_path):
    _, cfg = make_cfg(libtype=0, paired=True)
    d = lib_dirs(cfg)
    log = str(tmp_path / "pipe.log")
    attila.step_join(cfg, log)
    attila.step_filter(cfg, log)
    attila.step_translate(cfg, log)
    assert (d["final"] / "vh_final.fqjoinaa.fasta").exists()


@pytest.mark.xfail(raises=AssertionError, reason=(
    "paired-end: step_join checks `'input1' in f` on the file path instead of the config key "
    "(attila.py L145), so FastQC reports for R0 go to FinalRound"))
def test_paired_raw_fastqc_goes_to_matching_round(make_cfg, fake_tools, tmp_path):
    _, cfg = make_cfg(libtype=0, paired=True)
    d = lib_dirs(cfg)
    attila.step_join(cfg, str(tmp_path / "pipe.log"))
    fastqc = [c for c in fake_tools.commands("fastqc") if c.startswith("fastqc")]
    assert fastqc[0] == f"fastqc -q -o {d['init']} {cfg['input1r1dir']}"
    assert fastqc[1] == f"fastqc -q -o {d['init']} {cfg['input1r2dir']}"
    assert fastqc[2] == f"fastqc -q -o {d['final']} {cfg['input2r1dir']}"


@pytest.mark.xfail(raises=AssertionError, reason=(
    "paired-end: output_csv_counting counts 'X.fqjoin.fq' (attila.py L1178) but fastq-join writes 'X.fqjoin'"))
def test_paired_end_joining_count(make_cfg, fake_tools, tmp_path):
    _, cfg = make_cfg(libtype=0, paired=True)
    d = lib_dirs(cfg)
    attila.step_join(cfg, str(tmp_path / "pipe.log"))
    attila.output_csv_counting(cfg)
    rows = (d["lib"] / "vhSequenceCounting.csv").read_text().splitlines()
    assert "R0,12,joining" in rows
    assert "RN,28,joining" in rows


@pytest.mark.xfail(raises=ZeroDivisionError, reason=(
    "calculate_z_test divides by the library sizes without checking for 0 (attila.py L837/L845)"))
def test_z_test_with_empty_initial_library(tmp_path):
    enriched = write(tmp_path / "list.fasta", ">c1|FRAME:1+|FOLD-CHANGE:2.0|P1:0.0|P2:35714285714285.72\nEVQL\n")
    numbered = write(tmp_path / "numbered.fasta", ">c1|FRAME:1+|FOLD-CHANGE:2.0\nEVQL\n")
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), 0, 28, 1, str(out))
    assert "p-value:" in out.read_text()


@pytest.mark.xfail(raises=ZeroDivisionError, reason=(
    "find_duplicates divides by the R0 frequency, which is 0.0 when the R0 library size is 0 (attila.py L474)"))
def test_find_duplicates_with_zero_initial_frequency(tmp_path):
    r0 = write(tmp_path / "r0.txt", "#a|FRAME:1+|[1-3]|0|0.0\n*XXX\n>a|FRAME:1+\nMXXX\n")
    rn = write(tmp_path / "rn.txt", "#b|FRAME:1+|[1-3]|20|50000000000000.0\n*XXX\n>b|FRAME:1+\nMXXX\n")
    out = tmp_path / "enriched.fasta"
    attila.find_duplicates(str(r0), str(rn), str(out))
    assert out.exists()


@pytest.mark.xfail(raises=AssertionError, reason=(
    "run.sh compile_report writes index.html but announces Report.html (run.sh L337 vs L348)"))
def test_run_sh_announces_the_generated_report():
    text = (REPO_ROOT / "run.sh").read_text(encoding="utf-8")
    success = next(line for line in text.splitlines() if "Report generated successfully" in line)
    assert "index.html" in success
