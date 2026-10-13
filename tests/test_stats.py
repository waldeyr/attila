"""One-tailed z-test for the difference of proportions (calculate_z_test)."""

import math

import pytest

import attila
from helpers import write

SCALE = 1e14


def _inputs(tmp_path, clones):
    """clones: list of (id, count_r0, size_r0, count_rn, size_rn)."""
    enriched, numbered = [], []
    for cid, c0, n0, cn, nn in clones:
        p1, p2 = c0 / n0 * SCALE, cn / nn * SCALE
        enriched.append(f">{cid}|FRAME:1+|FOLD-CHANGE:{p2 / p1:.4f}|P1:{p1}|P2:{p2}\nEVQL\n")
        numbered.append(f">{cid}|FRAME:1+|FOLD-CHANGE:{p2 / p1:.4f}\nEVQL\n")
    return (write(tmp_path / "list.fasta", "".join(enriched)),
            write(tmp_path / "numbered.fasta", "".join(numbered)))


def parse_stats(path):
    lines = path.read_text().splitlines()
    alfa = float(lines[0].split("\t")[1])
    out, cur = {}, None
    for line in lines[1:]:
        if line.startswith(">"):
            cur = line[1:]
            out[cur] = {}
        else:
            key, value = line.split(":\t")
            out[cur][key] = float(value)
    return alfa, out


def reference_z_test(c0, n0, cn, nn):
    """Textbook pooled two-proportion z-test, one-tailed (H1: p_after > p_before), plus 95% CI of p0 - pn."""
    p0, pn = c0 / n0, cn / nn
    pooled = (c0 + cn) / (n0 + nn)
    se = math.sqrt(pooled * (1 - pooled) * (1 / n0 + 1 / nn))
    z = (pn - p0) / se
    p_value = 0.5 * math.erfc(z / math.sqrt(2))
    se_ci = math.sqrt(p0 * (1 - p0) / n0 + pn * (1 - pn) / nn)
    return p_value, (p0 - pn) - 1.96 * se_ci, (p0 - pn) + 1.96 * se_ci


def test_matches_reference_example_run(tmp_path):
    # Clone A of the synthetic VH set: 1/12 reads before, 10/28 after.
    enriched, numbered = _inputs(tmp_path, [("read_VH_A_1", 1, 12, 10, 28)])
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), 12, 28, 3, str(out))
    alfa, stats = parse_stats(out)
    assert alfa == pytest.approx(0.05 / 3)
    assert stats["read_VH_A_1"]["p-value"] == pytest.approx(0.03776210445990613)
    assert stats["read_VH_A_1"]["infIC"] == pytest.approx(-0.5103568333051713)
    assert stats["read_VH_A_1"]["supIC"] == pytest.approx(-0.03726221431387644)


@pytest.mark.parametrize("c0, n0, cn, nn", [
    (10, 10000, 200, 10000),
    (5, 1000, 6, 2000),
    (1, 12, 6, 28),
    (300, 50000, 900, 40000),
])
def test_matches_textbook_formula(tmp_path, c0, n0, cn, nn):
    enriched, numbered = _inputs(tmp_path, [("clone", c0, n0, cn, nn)])
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), n0, nn, 1, str(out))
    _, stats = parse_stats(out)
    p_value, inf_ic, sup_ic = reference_z_test(c0, n0, cn, nn)
    assert stats["clone"]["p-value"] == pytest.approx(p_value, rel=1e-9, abs=1e-12)
    assert stats["clone"]["infIC"] == pytest.approx(inf_ic, rel=1e-9, abs=1e-12)
    assert stats["clone"]["supIC"] == pytest.approx(sup_ic, rel=1e-9, abs=1e-12)


def test_strong_enrichment_is_significant_after_bonferroni(tmp_path):
    enriched, numbered = _inputs(tmp_path, [("hit", 10, 10000, 200, 10000), ("noise", 50, 10000, 52, 10000)])
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), 10000, 10000, 2, str(out))
    alfa, stats = parse_stats(out)
    assert alfa == pytest.approx(0.025)
    assert stats["hit"]["p-value"] < alfa < stats["noise"]["p-value"]
    assert stats["hit"]["supIC"] < 0 < stats["noise"]["supIC"]


def test_output_follows_numbered_order(tmp_path):
    enriched, numbered = _inputs(tmp_path, [("c1", 1, 100, 5, 100), ("c2", 1, 100, 9, 100), ("c3", 1, 100, 2, 100)])
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), 100, 100, 3, str(out))
    assert [l for l in out.read_text().splitlines() if l.startswith(">")] == [">c1", ">c2", ">c3"]


def test_degenerate_pooled_proportion_gives_p_value_one(tmp_path):
    # Both proportions are 1, so the pooled proportion is 1 and the standard error is 0.
    enriched, numbered = _inputs(tmp_path, [("all", 10, 10, 20, 20)])
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), 10, 20, 1, str(out))
    assert parse_stats(out)[1]["all"]["p-value"] == 1.0


def test_numbered_clone_missing_from_enriched_list(tmp_path):
    enriched = write(tmp_path / "list.fasta", "")
    numbered = write(tmp_path / "numbered.fasta", ">ghost|FRAME:1+|FOLD-CHANGE:2.0\nEVQL\n")
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), 10, 10, 1, str(out))
    assert parse_stats(out)[1]["ghost"] == {"p-value": 1.0, "infIC": 0.0, "supIC": 0.0}


def test_no_candidates(tmp_path):
    enriched = write(tmp_path / "list.fasta", "")
    numbered = write(tmp_path / "numbered.fasta", "")
    out = tmp_path / "stats.txt"
    attila.calculate_z_test(str(enriched), str(numbered), 10, 10, 0, str(out))
    assert out.read_text() == "No candidates\n"


def test_count_fasta_records(tmp_path):
    assert attila.count_fasta_records(str(write(tmp_path / "a.fasta", ">1\nA\n>2\nC\nG\n>3\n"))) == 3
    assert attila.count_fasta_records(str(tmp_path / "missing.fasta")) == 0
