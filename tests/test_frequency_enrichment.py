"""translate_all -> frequency_counter -> find_duplicates -> get_nsequences."""

import re

import pytest

import attila
from helpers import GENCODE, RN_COUNTS, SYNTHETIC_DIR, VARIANT_RESIDUE, read_fasta, synthetic_read, write

SCALE = 1e14  # frequencies are stored as count / libsize * 1e14
ENRICHED_HEADER = re.compile(r"^(.+)\|FOLD-CHANGE:([\d.]+)\|P1:([\d.e+-]+)\|P2:([\d.e+-]+)$")


def parse_enriched(path):
    out = []
    for header, seq in read_fasta(path):
        m = ENRICHED_HEADER.match(header)
        assert m, header
        out.append((m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4)), seq))
    return out


# --- translate_all ------------------------------------------------------------

def test_translate_all_output_format(tmp_path, codon_table):
    a = synthetic_read("vh_final", "read_VH_A_1")
    d = synthetic_read("vh_final", "read_VH_D_1")
    src = write(tmp_path / "in.fasta", f">read_A\n{a[:100]}\n{a[100:]}\n>junk\n{'ACGT' * 30}\n>read_D\n{d}\n")
    aa_out, nt_out = tmp_path / "aa.fasta", tmp_path / "nt.fasta"

    attila.translate_all(str(src), str(GENCODE), str(aa_out), str(nt_out), 0)

    prot_a = attila.translate_seq(a, 1, "+", codon_table)
    prot_d = attila.translate_seq(d, 1, "+", codon_table)
    assert aa_out.read_text().splitlines() == [
        ">read_A|FRAME:1+", prot_a, "#read_A|FRAME:1+|[97-105]", "AKGGGGGGG*",
        ">read_D|FRAME:1+", prot_d, "#read_D|FRAME:1+|[97-105]", "AKQQQQQQQ*",
        "<",
    ]
    # Reads without a recognisable V domain are dropped; multi-line FASTA is joined.
    assert read_fasta(nt_out) == [("read_A|FRAME:1+", a), ("read_D|FRAME:1+", d)]


def test_translate_all_empty_input_still_writes_end_marker(tmp_path):
    src = write(tmp_path / "in.fasta", "")
    aa_out = tmp_path / "aa.fasta"
    attila.translate_all(str(src), str(GENCODE), str(aa_out), str(tmp_path / "nt.fasta"), 1)
    assert aa_out.read_text() == "<\n"


# --- frequency_counter --------------------------------------------------------

def test_frequency_counter_on_synthetic_final_round(tmp_path):
    aa = tmp_path / "aa.fasta"
    filtered = SYNTHETIC_DIR / "vh_finalfiltered.fasta"
    attila.translate_all(str(filtered), str(GENCODE), str(aa), str(tmp_path / "nt.fasta"), 0)
    freq = tmp_path / "freq.txt"

    attila.frequency_counter(str(aa), str(filtered), str(freq))

    headers = [l for l in freq.read_text().splitlines() if l.startswith("#")]
    cdrs = [l[1:] for l in freq.read_text().splitlines() if l.startswith("*")]
    # One block per CDR3, most frequent first: A(10) D(8) E(6) C(3) B(1).
    order = sorted(RN_COUNTS, key=RN_COUNTS.get, reverse=True)
    assert cdrs == ["AK" + VARIANT_RESIDUE[v] * 7 for v in order]
    for header, variant in zip(headers, order):
        idcdr, libsize, value = header[1:].rsplit("|", 2)
        assert idcdr == f"read_VH_{variant}_1|FRAME:1+|[97-105]"
        assert int(libsize) == 28
        assert float(value) == pytest.approx(RN_COUNTS[variant] / 28 * SCALE)
    # Every read of a CDR3 is listed under it.
    block_a = freq.read_text().split("*AKGGGGGGG\n")[1].split("#")[0].splitlines()
    assert [l for l in block_a if l.startswith(">")] == [f">read_VH_A_{i}|FRAME:1+" for i in range(1, 11)]


def test_frequency_counter_without_filtered_file_reports_zero(tmp_path):
    aa = write(tmp_path / "aa.fasta", ">r1|FRAME:1+\nMAAA\n#r1|FRAME:1+|[2-4]\nAAA*\n<\n")
    freq = tmp_path / "freq.txt"
    attila.frequency_counter(str(aa), str(tmp_path / "missing.fasta"), str(freq))
    assert freq.read_text() == "#r1|FRAME:1+|[2-4]|0|0.0\n*AAA\n>r1|FRAME:1+\nMAAA\n"


def test_frequency_counter_stops_at_end_marker(tmp_path):
    aa = write(tmp_path / "aa.fasta",
               ">r1|FRAME:1+\nMAAA\n#r1|FRAME:1+|[2-4]\nAAA*\n<\n"
               ">r2|FRAME:1+\nMCCC\n#r2|FRAME:1+|[2-4]\nCCC*\n")
    filtered = write(tmp_path / "f.fasta", ">r1\nA\n>r2\nC\n")
    freq = tmp_path / "freq.txt"
    attila.frequency_counter(str(aa), str(filtered), str(freq))
    assert [l for l in freq.read_text().splitlines() if l.startswith("*")] == ["*AAA"]


# --- find_duplicates ----------------------------------------------------------

def _freq_block(idcdr, libsize, freq, cdr, reads):
    lines = [f"#{idcdr}|FRAME:1+|[1-3]|{libsize}|{freq}", f"*{cdr}"]
    for rid, seq in reads:
        lines += [f">{rid}|FRAME:1+", seq]
    return "\n".join(lines) + "\n"


@pytest.fixture
def freq_files(tmp_path):
    r0 = write(tmp_path / "r0freq.txt",
               _freq_block("a1", 10, 50000000000000.0, "YYY", [("a1", "MYYY")]) +
               _freq_block("a2", 10, 20000000000000.0, "XXX", [("a2", "MXXX")]))
    rn = write(tmp_path / "rnfreq.txt",
               _freq_block("b1", 20, 80000000000000.0, "XXX", [("b1", "MAXXX"), ("b2", "MAAXXX"), ("b3", "MXXX")]) +
               _freq_block("b4", 20, 50000000000000.0, "WWW", [("b4", "MWWW")]) +
               _freq_block("b5", 20, 15000000000000.0, "ZZZ", [("b5", "MZZZ")]) +
               _freq_block("b6", 20, 10000000000000.0, "YYY", [("b6", "MYYY")]) +
               _freq_block("b7", 20, 10000000000000.0, "VVV", [("b7", "MVVV")]))
    return r0, rn


def test_find_duplicates_fold_change_and_ranking(tmp_path, freq_files):
    out = tmp_path / "enriched.fasta"
    attila.find_duplicates(str(freq_files[0]), str(freq_files[1]), str(out))
    got = parse_enriched(out)

    # WWW and ZZZ are absent from R0, so their R0 frequency is 1/libsize(R0) = 1e13.
    # YYY decreased and VVV is unchanged at the pseudo-count, so both are dropped.
    assert [(rid, fold) for rid, fold, *_ in got] == [
        ("b4|FRAME:1+", 5.0),
        ("b2|FRAME:1+", 4.0),
        ("b5|FRAME:1+", 1.5),
    ]
    assert got[1][2:4] == pytest.approx((2e13, 8e13))
    assert got[0][2] == pytest.approx(1e13)
    # For a CDR3 with several reads, the longest sequence represents it.
    assert got[1][4] == "MAAXXX"


def test_find_duplicates_fold_change_has_four_decimals(tmp_path):
    r0 = write(tmp_path / "r0.txt", _freq_block("a", 3, 30000000000000.0, "XXX", [("a", "MXXX")]))
    rn = write(tmp_path / "rn.txt", _freq_block("b", 7, 100000000000000.0, "XXX", [("b", "MXXX")]))
    out = tmp_path / "e.fasta"
    attila.find_duplicates(str(r0), str(rn), str(out))
    assert "|FOLD-CHANGE:3.3333|" in out.read_text()


# --- get_nsequences -----------------------------------------------------------

@pytest.mark.parametrize("n, expected", [(1, ["r1"]), (2, ["r1", "r2"]), (3, ["r1", "r2", "r3"]), (10, ["r1", "r2", "r3"])])
def test_get_nsequences(tmp_path, n, expected):
    src = write(tmp_path / "in.fasta", ">r1\nAA\nAA\n>r2\nCC\n>r3\nGG\n")
    out = tmp_path / "out.fasta"
    attila.get_nsequences(str(src), str(out), n)
    records = read_fasta(out)
    assert [h for h, _ in records] == expected
    assert records[0][1] == "AAAA"
