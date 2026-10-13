"""Read counting, the SequenceCounting.csv summary, and .cfg parsing."""

import pytest

import attila
from helpers import SYNTHETIC_DIR, lib_dirs, write


# --- counting -----------------------------------------------------------------

def test_count_sequences_step_fastq_counts_records():
    assert attila.count_sequences_step(str(SYNTHETIC_DIR / "vh_initial.fastq"), is_fastq=True) == 12
    assert attila.count_sequences_step(str(SYNTHETIC_DIR / "vh_final.fastq"), is_fastq=True) == 28


def test_count_sequences_step_fasta_counts_headers(tmp_path):
    fasta = write(tmp_path / "x.fasta", ">a\nAC\nGT\n>b\nAA\n")
    assert attila.count_sequences_step(str(fasta), is_fastq=False) == 2


@pytest.mark.parametrize("is_fastq", [True, False])
def test_count_sequences_step_missing_file(tmp_path, is_fastq):
    assert attila.count_sequences_step(str(tmp_path / "missing"), is_fastq=is_fastq) == 0


def test_count_cdr_headers(tmp_path):
    freq = write(tmp_path / "freq.txt", "#a|1|2\n*AAA\n>a\nM\n#b|1|2\n*CCC\n>b\nM\n")
    assert attila.count_cdr_headers(str(freq)) == 2
    assert attila.count_cdr_headers(str(tmp_path / "missing.txt")) == 0


def test_output_csv_counting_single_end(make_cfg):
    _, cfg = make_cfg(libtype=0)
    d = lib_dirs(cfg)
    write(d["init"] / "vh_initialfiltered.fasta", ">a\nA\n>b\nA\n>c\nA\n")
    write(d["final"] / "vh_finalfiltered.fasta", ">a\nA\n>b\nA\n")
    write(d["init"] / "vh_initialaa.fasta", ">a\nM\n#a\nX*\n>b\nM\n#b\nX*\n<\n")
    write(d["final"] / "vh_finalaa.fasta", ">a\nM\n#a\nX*\n<\n")
    write(d["init"] / "vh_initialfreq.txt", "#a|3|1\n*X\n")
    write(d["final"] / "vh_finalfreq.txt", "#a|2|1\n*X\n#b|2|1\n*Y\n")
    write(d["selected"] / "vhenriched.fasta", ">e1\nM\n>e2\nM\n>e3\nM\n>e4\nM\n")
    write(d["selected"] / "vhlist5numbered.fasta", ">e1\nM\n")

    attila.output_csv_counting(cfg)

    assert (d["lib"] / "vhSequenceCounting.csv").read_text().splitlines() == [
        "library,reads,step",
        "R0,12,raw", "RN,28,raw",
        "R0,3,filtering", "RN,2,filtering",
        "R0,2,translation", "RN,1,translation",
        "R0,1,frequency", "RN,2,frequency",
        "Selected,4,enrichment",
        "Selected,1,numeration",
    ]


def test_output_csv_counting_vl_with_nothing_processed(make_cfg):
    _, cfg = make_cfg(libtype=1, listsize=10)
    d = lib_dirs(cfg)
    attila.output_csv_counting(cfg)
    rows = (d["lib"] / "vlSequenceCounting.csv").read_text().splitlines()
    assert rows[1:3] == ["R0,12,raw", "RN,28,raw"]
    assert all(row.split(",")[1] == "0" for row in rows[3:])
    assert "joining" not in "".join(rows)


# --- parse_cfg ----------------------------------------------------------------

def test_parse_cfg(tmp_path):
    cfg_file = write(tmp_path / "x.cfg",
                     "# comment: ignored\n"
                     "----------------\n"
                     "\n"
                     "projectname:   my project  \n"
                     "projectdir: /data/runs\n"
                     "igblastdir: C:/tools/igblast\n"
                     "libtype:0\n"
                     "no separator here\n")
    assert attila.parse_cfg(str(cfg_file)) == {
        "projectname": "my project",
        "projectdir": "/data/runs",
        "igblastdir": "C:/tools/igblast",
        "libtype": "0",
    }


def test_parse_cfg_shipped_templates():
    cfg = attila.parse_cfg(str(SYNTHETIC_DIR / "example_VL.cfg"))
    assert cfg["libtype"] == "1"
    assert cfg["listsize"] == "5"
    assert cfg["input1dir"].endswith("vl_initial.fastq")
    assert cfg["projectdir"] == "@WORKSPACE_DIR@/example_output"


def test_parse_cfg_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        attila.parse_cfg(str(tmp_path / "missing.cfg"))
