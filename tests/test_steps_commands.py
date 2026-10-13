"""Steps that shell out to external tools, checked through the commands they build (fake_tools)."""

import pytest

import attila
from helpers import SYNTHETIC_DIR, lib_dirs, read_fasta, write


def test_step_join_single_end_runs_fastqc_per_round(make_cfg, fake_tools, tmp_path):
    _, cfg = make_cfg(libtype=0)
    d = lib_dirs(cfg)
    attila.step_join(cfg, str(tmp_path / "pipe.log"))

    assert fake_tools.commands("fastqc") == [
        f"fastqc -q -o {d['init']} {cfg['input1dir']}",
        f"fastqc -q -o {d['final']} {cfg['input2dir']}",
    ]
    assert not fake_tools.commands("fastq-join")


def test_step_join_paired_end_joins_both_rounds(make_cfg, fake_tools, tmp_path):
    _, cfg = make_cfg(libtype=1, paired=True)
    d = lib_dirs(cfg)
    attila.step_join(cfg, str(tmp_path / "pipe.log"))

    joins = [c for c in fake_tools.commands("-o ") if not c.startswith("fastqc")]
    assert len(joins) == 2
    assert f"{cfg['input1r1dir']} {cfg['input1r2dir']} -o {d['init']}/vl_initial.fq" in joins[0]
    assert f"{cfg['input2r1dir']} {cfg['input2r2dir']} -o {d['final']}/vl_final.fq" in joins[1]
    assert len([c for c in fake_tools.commands("fastqc ") if c.startswith("fastqc")]) == 4
    assert (d["init"] / "vl_initial.fqjoin").exists()


def test_step_filter_builds_prinseq_commands(make_cfg, fake_tools, tmp_path):
    _, cfg = make_cfg(libtype=0, minlen=320, minqual=25)
    d = lib_dirs(cfg)
    attila.step_filter(cfg, str(tmp_path / "pipe.log"))

    prinseq = fake_tools.commands("prinseq-lite")
    assert len(prinseq) == 4
    before, filtered = prinseq[:2], prinseq[2:]
    assert all("-out_format 1" in c and "-min_len" not in c for c in before)
    assert all("-out_format 5" in c and "-min_qual_mean 25" in c and "-min_len 320" in c for c in filtered)
    assert f"-out_good {d['init']}/vh_initialfiltered " in filtered[0]
    assert f"-out_good {d['final']}/vh_finalfiltered " in filtered[1]
    assert len(read_fasta(d["final"] / "vh_finalfiltered.fasta")) == 28
    assert fake_tools.commands("fastqc -q -o")[-1].endswith(f"{d['final']}/vh_finalfiltered.fastq >> {tmp_path}/pipe.log 2>&1")


def test_step_filter_without_thresholds(make_cfg, fake_tools, tmp_path):
    _, cfg = make_cfg(libtype=0, minlen="", minqual="")
    lib_dirs(cfg)
    attila.step_filter(cfg, str(tmp_path / "pipe.log"))
    assert not any("-min_" in c for c in fake_tools.commands("prinseq-lite"))


def test_step_filter_sanitizes_illumina_headers(make_cfg, fake_tools, tmp_path):
    seq = "ACGT" * 80
    qual = "I" * len(seq)
    raw = write(tmp_path / "raw" / "r0.fastq",
                f"@M00123:45:000000000-ABCDE:1:1101:15589:1331 1:N:0:1\n{seq}\n+\n{qual}\n"
                f"@plainread\n{seq}\n+\n{qual}\n")
    _, cfg = make_cfg(libtype=0, input1dir=raw)
    d = lib_dirs(cfg)

    attila.step_filter(cfg, str(tmp_path / "pipe.log"))

    sanitized = d["init"] / "r0newid.fq"
    assert sanitized.read_text().splitlines()[0] == "@M00123:45:000000000-ABCDE:1:1101:15589:1331"
    assert sanitized.read_text().splitlines()[4] == "@plainread"
    assert all(f"-fastq {sanitized} " in c for c in fake_tools.commands("prinseq-lite")[::2])
    # The final-round file has no spaces in its headers, so it is used as is.
    assert not (d["final"] / "vh_finalnewid.fq").exists()


@pytest.mark.parametrize("libtype, lib", [(0, "VH"), (1, "VL")])
def test_step_germline_command(make_cfg, fake_tools, tmp_path, monkeypatch, libtype, lib):
    monkeypatch.setattr(attila.shutil, "which", lambda name: "/opt/ncbi-igblast-1.22.0/bin/igblastp")
    _, cfg = make_cfg(libtype=libtype)
    d = lib_dirs(cfg)
    numbered = write(d["selected"] / f"{lib.lower()}list5numbered.fasta", ">c1|FRAME:1+|FOLD-CHANGE:2.0\nEVQL\n")

    attila.step_germline(cfg, str(tmp_path / "pipe.log"))

    cmd, kwargs = fake_tools.calls[-1]
    assert cmd.startswith("/opt/ncbi-igblast-1.22.0/bin/igblastp ")
    assert f"-germline_db_V {cfg['packagedir']}/data/germline_human/{lib}_germline_human " in cmd
    assert f"-query {numbered} " in cmd
    assert "-organism human -domain_system kabat -outfmt 7" in cmd
    assert cmd.endswith(f"> {d['selected']}/{lib.lower()}list5numberedgermlineclassification.txt")
    assert kwargs["check"] is True
    assert kwargs["env"]["IGDATA"] == "/opt/ncbi-igblast-1.22.0"


def test_step_germline_falls_back_to_igblastdir(make_cfg, fake_tools, tmp_path, monkeypatch):
    monkeypatch.setattr(attila.shutil, "which", lambda name: None)
    _, cfg = make_cfg(libtype=0, igblastdir=tmp_path / "igblast")
    d = lib_dirs(cfg)
    write(d["selected"] / "vhlist5numbered.fasta", ">c1|FRAME:1+|FOLD-CHANGE:2.0\nEVQL\n")
    attila.step_germline(cfg, str(tmp_path / "pipe.log"))
    cmd, kwargs = fake_tools.calls[-1]
    assert cmd.startswith(f"{tmp_path}/igblast/bin/igblastp ")
    assert kwargs["env"]["IGDATA"] == str(tmp_path / "igblast")


# --- report step: plotting -----------------------------------------------------

def _report_inputs(cfg):
    d = lib_dirs(cfg)
    write(d["init"] / "vh_initialbeforefiltering.fasta", ">a\n" + "A" * 300 + "\n>b\n" + "A" * 10 + "\n")
    write(d["final"] / "vh_finalbeforefiltering.fasta", ">a\n" + "A" * 300 + "\n")
    write(d["lib"] / "vhSequenceCounting.csv",
          "library,reads,step\nR0,2,raw\nRN,1,raw\nR0,1,filtering\nRN,1,filtering\nSelected,1,enrichment\n")
    return d


def test_step_report_r_fallback_without_matplotlib(make_cfg, fake_tools, tmp_path, monkeypatch):
    monkeypatch.setitem(__import__("sys").modules, "matplotlib", None)  # import matplotlib -> ImportError
    _, cfg = make_cfg(libtype=0)
    d = _report_inputs(cfg)

    attila.step_report(cfg, str(tmp_path / "pipe.log"))

    plot1 = (d["report"] / "plot1vh.r").read_text()
    assert "Initial\t\t\t\t\tAdequate\t\t\t\t\t0.5000" in plot1
    assert "Final\t\t\t\t\tAdequate\t\t\t\t\t1.0000" in plot1
    assert f'png(file="{d["report"]}/plot1vh.png"' in plot1
    plot2 = (d["report"] / "plot2vh.r").read_text()
    assert "Initial\t\t\t\t\tNone\t\t\t\t\t2" in plot2
    assert "Final\t\t\t\t\tEnrichment\t\t\t\t\t1" in plot2
    assert fake_tools.commands("Rscript") == [f"Rscript {d['report']}/plot1vh.r", f"Rscript {d['report']}/plot2vh.r"]
    # The unfiltered FASTA files are deleted once plotted.
    assert not (d["init"] / "vh_initialbeforefiltering.fasta").exists()
    assert not (d["final"] / "vh_finalbeforefiltering.fasta").exists()


def test_step_report_with_matplotlib(make_cfg, fake_tools, tmp_path):
    pytest.importorskip("matplotlib")
    pytest.importorskip("numpy")
    _, cfg = make_cfg(libtype=0)
    d = _report_inputs(cfg)

    attila.step_report(cfg, str(tmp_path / "pipe.log"))

    for name in ("plot1vh.png", "plot2vh.png"):
        assert (d["report"] / name).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert not fake_tools.commands("Rscript")
