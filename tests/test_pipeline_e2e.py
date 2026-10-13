"""End-to-end runs of attila.main() + html_creator.main() on the synthetic dataset.

External tools and the Abnum web service are replaced by fake_tools / fake_abnum,
so this checks the pipeline's own logic. The expected values follow from the
dataset's design: clones A, D and E are enriched from R0 to RN, B decreases and
C stays the same.
"""

import shutil
import sys

import pytest

import attila
import html_creator
from helpers import SYNTHETIC_DIR, read_fasta, synthetic_read

pytestmark = pytest.mark.integration

LIBS = [(0, "VH"), (1, "VL")]


def run_attila(monkeypatch, cfg_path, tmp_path, steps=None):
    if steps:
        monkeypatch.setenv("ATTILA_STEPS", steps)
    monkeypatch.setattr(sys, "argv", ["attila.py", str(cfg_path), str(tmp_path / "attila.log")])
    attila.main()


@pytest.mark.parametrize("libtype, lib", LIBS)
def test_full_pipeline(make_cfg, fake_tools, fake_abnum, tmp_path, monkeypatch, libtype, lib):
    cfg_path, cfg = make_cfg(libtype=libtype)
    aux = lib.lower()

    run_attila(monkeypatch, cfg_path, tmp_path)

    base = tmp_path / "out" / "proj" / lib
    selected = base / "SelectedSequences"

    # Enrichment: exactly A, D, E, ranked by fold change (RN/R0 frequency ratios).
    enriched = read_fasta(selected / f"{aux}enriched.fasta")
    assert [h.split("|P1")[0] for h, _ in enriched] == [
        f"read_{lib}_A_1|FRAME:1+|FOLD-CHANGE:4.2857",   # (10/28) / (1/12)
        f"read_{lib}_D_1|FRAME:1+|FOLD-CHANGE:3.4286",   # (8/28)  / (1/12)
        f"read_{lib}_E_1|FRAME:1+|FOLD-CHANGE:1.2857",   # (6/28)  / (2/12)
    ]
    assert len(read_fasta(selected / f"{aux}list5.fasta")) == 3

    # Numbering went through Abnum once per candidate and kept the protein intact.
    assert len(fake_abnum.queried) == 3
    numbered = read_fasta(selected / f"{aux}list5numbered.fasta")
    assert [seq for _, seq in numbered] == [seq for _, seq in enriched]

    # Nucleotide recovery returns the original reads.
    nt = read_fasta(selected / f"{aux}list5numberednt.fasta")
    assert nt[0] == (f"read_{lib}_A_1|FOLD-CHANGE:4.2857", synthetic_read(f"{aux}_final", f"read_{lib}_A_1"))

    # Statistics: Bonferroni over 3 candidates; clone A has the smallest p-value.
    stats = (tmp_path / "out" / "proj" / "Report" / f"{aux}outputRstats.txt").read_text().splitlines()
    assert stats[0] == f"alfaBonf:\t{0.05 / 3}"
    assert stats[1] == f">read_{lib}_A_1"
    assert float(stats[2].split("\t")[1]) == pytest.approx(0.03776210445990613)

    # Read accounting.
    assert (base / f"{aux}SequenceCounting.csv").read_text().splitlines() == [
        "library,reads,step",
        "R0,12,raw", "RN,28,raw",
        "R0,12,filtering", "RN,28,filtering",
        "R0,12,translation", "RN,28,translation",
        "R0,5,frequency", "RN,5,frequency",
        "Selected,3,enrichment",
        "Selected,3,numeration",
    ]

    # External tools were called for every step that needs them.
    assert fake_tools.commands("prinseq-lite")
    assert fake_tools.commands("igblastp -germline_db_V")
    assert (selected / f"{aux}list5numberedgermlineclassification.txt").exists()
    # The report step removes the unfiltered FASTA files.
    assert not list(base.rglob("*beforefiltering.fasta"))


def test_steps_env_selects_steps(make_cfg, fake_tools, fake_abnum, tmp_path, monkeypatch):
    cfg_path, cfg = make_cfg(libtype=0)
    init = tmp_path / "out" / "proj" / "VH" / "InitialRound"
    final = tmp_path / "out" / "proj" / "VH" / "FinalRound"
    init.mkdir(parents=True)
    final.mkdir(parents=True)
    shutil.copy(SYNTHETIC_DIR / "vh_initialfiltered.fasta", init)
    shutil.copy(SYNTHETIC_DIR / "vh_finalfiltered.fasta", final)

    run_attila(monkeypatch, cfg_path, tmp_path, steps="Translate, frequency")

    selected = tmp_path / "out" / "proj" / "VH" / "SelectedSequences"
    assert fake_tools.calls == []
    assert fake_abnum.queried == []
    assert (final / "vh_finalaa.fasta").exists()
    assert (selected / "vhenriched.fasta").exists()
    assert not (selected / "vhlist5numbered.fasta").exists()
    # The counting CSV is always written.
    assert (tmp_path / "out" / "proj" / "VH" / "vhSequenceCounting.csv").exists()


def test_main_creates_project_layout(make_cfg, fake_tools, tmp_path, monkeypatch):
    cfg_path, _ = make_cfg(libtype=1)
    run_attila(monkeypatch, cfg_path, tmp_path, steps="join")
    proj = tmp_path / "out" / "proj"
    for d in ("VL/InitialRound", "VL/FinalRound", "VL/SelectedSequences", "Report"):
        assert (proj / d).is_dir()
    assert (proj / "vlPipelineError.log").exists()
    assert len(fake_tools.commands("fastqc")) == 2


def test_vh_and_vl_then_html_report(make_cfg, fake_tools, fake_abnum, tmp_path, monkeypatch):
    for libtype, _ in LIBS:
        cfg_path, _ = make_cfg(libtype=libtype)
        run_attila(monkeypatch, cfg_path, tmp_path)

    proj = tmp_path / "out" / "proj"
    sel = {lib: proj / lib / "SelectedSequences" for _, lib in LIBS}
    monkeypatch.setattr(sys, "argv", [
        "html_creator.py",
        str(sel["VH"] / "vhlist5numbered.fasta"), str(sel["VL"] / "vllist5numbered.fasta"),
        str(sel["VH"] / "vhlist5numberedgermlineclassification.txt"),
        str(sel["VL"] / "vllist5numberedgermlineclassification.txt"),
        str(proj / "index.html"),
        "Report/plot1vh.png", "Report/plot2vh.png", "Report/plot1vl.png", "Report/plot2vl.png",
        str(proj / "VH" / "vhSequenceCounting.csv"), str(proj / "VL" / "vlSequenceCounting.csv"),
        str(proj / "Report" / "vhoutputRstats.txt"), str(proj / "Report" / "vloutputRstats.txt"),
    ])
    html_creator.main()

    html = (proj / "index.html").read_text(encoding="utf-8")
    for lib in ("VH", "VL"):
        for clone in "ADE":
            assert f'data-bs-content="read_{lib}_{clone}_1"' in html
    assert "<td>4.29</td>" in html and "<td>3.43</td>" in html and "<td>1.29</td>" in html
    assert "<td>VH3-23</td>" in html and "<td>V1-13</td>" in html
    assert "No candidate clones found" not in html
    assert 'href="vhPipelineError.log"' in html
