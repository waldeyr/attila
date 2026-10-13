"""html_creator: input parsers and the generated index.html."""

import sys

import pytest

import html_creator
from helpers import write

VH_PROTEIN = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGGGGGGGWGQGTLVTVSS"


# --- parse_fasta --------------------------------------------------------------

def test_parse_fasta_continues_counter_across_libraries(tmp_path):
    vh = write(tmp_path / "vh.fasta", ">vh1|FRAME:1+|FOLD-CHANGE:4.2857\nEVQ-L\nAK\n>vh2|FRAME:2-|FOLD-CHANGE:1.5\nQVQ\n")
    vl = write(tmp_path / "vl.fasta", ">not a candidate header\nXXX\n>vl1|FRAME:1+|FOLD-CHANGE:10.0\nQSV\n")
    table1, table2 = {}, {}

    count_vh = html_creator.parse_fasta(str(vh), table1, table2, 0)
    count_vl = html_creator.parse_fasta(str(vl), table1, table2, count_vh)

    assert (count_vh, count_vl) == (2, 3)
    assert {k: (v["fc"], v["count"]) for k, v in table1.items()} == {
        "vh1": ("4.29", 0), "vh2": ("1.50", 1), "vl1": ("10.00", 2)}
    # Alignment gaps are removed and multi-line sequences joined.
    assert table2["vh1"]["seq"] == "EVQLAK"
    assert table2["vl1"]["seq"] == "QSV"


def test_parse_fasta_missing_file_keeps_counter(tmp_path):
    assert html_creator.parse_fasta(str(tmp_path / "missing"), {}, {}, 7) == 7


# --- parse_igblast ------------------------------------------------------------

def test_parse_igblast_regions_and_germline(tmp_path):
    numbered = write(tmp_path / "vh.fasta", f">c1|FRAME:1+|FOLD-CHANGE:4.2857\n{VH_PROTEIN}\n")
    igblast = write(tmp_path / "igblast.txt", "\n".join([
        "# IGBLASTP 1.22.0",
        "# Query: c1|FRAME:1+|FOLD-CHANGE:4.2857",
        "# Database: VH_germline_human",
        "FR1\t1\t25\t25\t25\t0\t0\t100",
        "CDR1\t26\t35\t10\t10\t0\t0\t100",
        "FR2\t36\t49\t14\t14\t0\t0\t100",
        "CDR2\t50\t66\t17\t17\t0\t0\t100",
        "FR3\t67\t98\t32\t32\t0\t0\t100",
        "V\tc1|FRAME:1+|FOLD-CHANGE:4.2857\tIGHV3-23*01\t98.2\t98\t2",
        "V\tc1|FRAME:1+|FOLD-CHANGE:4.2857\tIGHV3-23*04\t97.0\t98\t3",
        "# Query: unknown|FRAME:1+|FOLD-CHANGE:1.0",
        "FR1\t1\t25",
    ]) + "\n")
    table1, table2 = {}, {}
    html_creator.parse_fasta(str(numbered), table1, table2, 0)

    html_creator.parse_igblast(str(igblast), table1, table2)

    regions = table2["c1"]
    assert regions["fr1"] == VH_PROTEIN[0:25]
    assert regions["cdr1"] == "GFTFSSYAMS"
    assert regions["fr3"] == VH_PROTEIN[66:98]
    assert regions["cdr3"] == "GGGGGGG"  # Kabat FR3 ends at ...YYCAK
    assert regions["fr4"] == "WGQGTLVTVSS"
    # Only the best V hit (first line) is kept.
    assert (table1["c1"]["germline"], table1["c1"]["identity"]) == ("IGHV3-23*01", 98.2)
    assert "unknown" not in table2


def test_parse_igblast_missing_file(tmp_path):
    table1 = {"c1": {"germline": ""}}
    html_creator.parse_igblast(str(tmp_path / "missing"), table1, {})
    assert table1 == {"c1": {"germline": ""}}


# --- parse_csv_loss / parse_stats ---------------------------------------------

def test_parse_csv_loss(tmp_path):
    csv = write(tmp_path / "vh.csv",
                "library,reads,step\nR0,1000,raw\nRN,2000,raw\nR0,800,filtering\nRN,1500,filtering\nR0,700,translation\n")
    loss = {}
    html_creator.parse_csv_loss(str(csv), "vh", loss)
    html_creator.parse_csv_loss(str(tmp_path / "missing.csv"), "vl", loss)
    assert loss == {"vh": {"r0": {"raw": 1000, "filtering": 800}, "rn": {"raw": 2000, "filtering": 1500}}}


def test_parse_stats(tmp_path):
    stats_file = write(tmp_path / "stats.txt",
                       "alfaBonf:\t0.016666666666666666\n>c1\np-value:\t0.0377\ninfIC:\t-0.51\nsupIC:\t-0.037\n>c2\np-value:\t0.08\n")
    stats = {}
    assert html_creator.parse_stats(str(stats_file), stats) == "0.016666666666666666"
    assert stats == {"c1": {"pvalue": "0.0377", "infIC": "-0.51", "supIC": "-0.037"},
                     "c2": {"pvalue": "0.08", "infIC": "", "supIC": ""}}


def test_parse_stats_defaults(tmp_path):
    assert html_creator.parse_stats(str(tmp_path / "missing"), {}) == "0.05"
    assert html_creator.parse_stats(str(write(tmp_path / "s.txt", "No candidates\n")), {}) == "0.05"


# --- main ---------------------------------------------------------------------

def test_main_requires_13_arguments(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["html_creator.py", "a", "b"])
    with pytest.raises(SystemExit) as exc:
        html_creator.main()
    assert exc.value.code == 1
    assert "Usage:" in capsys.readouterr().out


def _run_main(monkeypatch, project, vh_numbered, vl_numbered):
    sel_vh, sel_vl, report = project / "VH" / "SelectedSequences", project / "VL" / "SelectedSequences", project / "Report"
    monkeypatch.setattr(sys, "argv", [
        "html_creator.py",
        str(vh_numbered), str(vl_numbered),
        str(sel_vh / "vhgerm.txt"), str(sel_vl / "vlgerm.txt"),
        str(project / "index.html"),
        "Report/plot1vh.png", "Report/plot2vh.png", "Report/plot1vl.png", "Report/plot2vl.png",
        str(project / "VH" / "vhSequenceCounting.csv"), str(project / "VL" / "vlSequenceCounting.csv"),
        str(report / "vhoutputRstats.txt"), str(report / "vloutputRstats.txt"),
    ])
    html_creator.main()
    return (project / "index.html").read_text(encoding="utf-8")


def test_main_writes_report(tmp_path, monkeypatch):
    project = tmp_path / "proj"
    vh = write(project / "VH" / "SelectedSequences" / "vhnum.fasta",
               f">vhA|FRAME:1+|FOLD-CHANGE:4.2857\n{VH_PROTEIN}\n>vhD|FRAME:1+|FOLD-CHANGE:3.4286\n{VH_PROTEIN}\n")
    vl = project / "VL" / "SelectedSequences" / "vlnum.fasta"  # no VL candidates
    write(project / "VH" / "SelectedSequences" / "vhgerm.txt",
          "# Query: vhA|FRAME:1+|FOLD-CHANGE:4.2857\nFR3\t67\t98\nV\tvhA|FRAME:1+|FOLD-CHANGE:4.2857\tIGHV3-23*01\t98.2\t1\n")
    write(project / "VH" / "vhSequenceCounting.csv", "library,reads,step\nR0,200,raw\nRN,400,raw\nR0,150,filtering\nRN,400,filtering\n")
    write(project / "Report" / "vhoutputRstats.txt", "alfaBonf:\t0.025\n>vhA\np-value:\t0.0012\ninfIC:\t-0.4\nsupIC:\t-0.1\n")
    write(project / "vherror.log", "")
    write(project / "Report" / "webpage.log", "")

    html = _run_main(monkeypatch, project, vh, vl)

    assert "<title>Selected sequences - ATTILA Report</title>" in html
    assert "VH Library" in html and "VL Library" in html
    # Loss of reads: (200-150)/200 and (400-400)/400; VL has no CSV.
    assert "<td>25.00%</td>\n<td>0.00%</td>" in html
    assert "<td>VL</td>\n<td>0.00%</td>\n<td>0.00%</td>" in html
    assert "&alpha;<sub>Bonf</sub>=0.025" in html
    assert 'data-bs-content="vhA">1\xb0</a></td>\n<td>4.29</td>\n<td>0.0012</td>\n<td>-0.4</td>\n<td>-0.1</td>\n<td>IGHV3-23*01</td>\n<td>98.2</td>' in html
    assert 'data-bs-content="vhD">2\xb0</a></td>\n<td>3.43</td>' in html
    assert '<td class="bg-cdr3">GGGGGGG</td>' in html
    assert '<td class="bg-fr">WGQGTLVTVSS</td>' in html
    assert html.count("No candidate clones found for this library.") == 1
    # Only logs that exist are linked.
    assert 'href="vherror.log"' in html and 'href="Report/webpage.log"' in html
    assert 'href="vlerror.log"' not in html


def test_main_without_any_candidates_or_logs(tmp_path, monkeypatch):
    project = tmp_path / "proj"
    html = _run_main(monkeypatch, project, project / "none_vh.fasta", project / "none_vl.fasta")
    assert html.count("No candidate clones found for this library.") == 2
    assert "No log files found." in html
    assert "&alpha;" not in html
