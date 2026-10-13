"""Opt-in checks against the real external tools: `pytest -m tools`.

Each test is skipped when its binary is not on PATH.
"""

import shutil

import pytest

import attila
import html_creator
from helpers import SYNTHETIC_DIR, lib_dirs, read_fasta, write

pytestmark = pytest.mark.tools


def require(*tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        pytest.skip(f"not on PATH: {', '.join(missing)}")


@pytest.mark.parametrize("libtype, prefix", [(0, "vh"), (1, "vl")])
def test_prinseq_filter_matches_reference_fasta(make_cfg, tmp_path, libtype, prefix):
    require("prinseq-lite", "fastqc")
    _, cfg = make_cfg(libtype=libtype)
    d = lib_dirs(cfg)

    attila.step_filter(cfg, str(tmp_path / "pipe.log"))

    for round_dir, name in ((d["init"], f"{prefix}_initial"), (d["final"], f"{prefix}_final")):
        assert read_fasta(round_dir / f"{name}filtered.fasta") == read_fasta(SYNTHETIC_DIR / f"{name}filtered.fasta")
        assert (round_dir / f"{name}beforefiltering.fasta").exists()


def test_prinseq_min_len_drops_short_reads(make_cfg, tmp_path):
    require("prinseq-lite", "fastqc")
    _, cfg = make_cfg(libtype=1, minlen=330)  # VL reads are 324 bp
    d = lib_dirs(cfg)
    attila.step_filter(cfg, str(tmp_path / "pipe.log"))
    # prinseq writes no output file at all when every read is filtered out.
    filtered = d["final"] / "vl_finalfiltered.fasta"
    assert not filtered.exists() or read_fasta(filtered) == []


def test_igblastp_classifies_synthetic_vh(make_cfg, tmp_path, codon_table):
    require("igblastp")
    _, cfg = make_cfg(libtype=0)
    d = lib_dirs(cfg)
    dna = read_fasta(SYNTHETIC_DIR / "vh_finalfiltered.fasta")[0][1]
    protein = attila.translate_seq(dna, 1, "+", codon_table)
    numbered = write(d["selected"] / "vhlist5numbered.fasta", f">read_VH_A_1|FRAME:1+|FOLD-CHANGE:4.2857\n{protein}\n")

    attila.step_germline(cfg, str(tmp_path / "pipe.log"))

    table1, table2 = {}, {}
    html_creator.parse_fasta(str(numbered), table1, table2, 0)
    html_creator.parse_igblast(str(d["selected"] / "vhlist5numberedgermlineclassification.txt"), table1, table2)
    # The synthetic framework is the human VH3-23 germline (named "VH3-23" in data/germline_human).
    assert table1["read_VH_A_1"]["germline"] == "VH3-23"
    assert table1["read_VH_A_1"]["identity"] == 100.0
    assert table2["read_VH_A_1"]["cdr1"]
