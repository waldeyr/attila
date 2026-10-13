"""Pure helpers: reverse complement, codon table, translation and CDR3 detection."""

import pytest

import attila
from helpers import GENCODE, VARIANT_RESIDUE, synthetic_read


# --- reverse_complement -------------------------------------------------------

@pytest.mark.parametrize("seq, expected", [
    ("ATCG", "CGAT"),
    ("atcg", "CGAT"),
    ("NnATCG", "CGATNN"),
    ("", ""),
])
def test_reverse_complement(seq, expected):
    assert attila.reverse_complement(seq) == expected


def test_reverse_complement_is_an_involution():
    seq = synthetic_read("vh_final", "read_VH_A_1")
    assert attila.reverse_complement(attila.reverse_complement(seq)) == seq


# --- codon table --------------------------------------------------------------

def test_load_codon_table_has_all_64_codons(codon_table):
    assert len(codon_table) == 64
    assert codon_table["ATG"] == "M"
    assert codon_table["TGG"] == "W"
    assert codon_table["TGT"] == codon_table["TGC"] == "C"


@pytest.mark.parametrize("stop", ["TAA", "TAG", "TGA"])
def test_stop_codons_translate_to_dash(codon_table, stop):
    assert codon_table[stop] == "-"


def test_load_codon_table_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        attila.load_codon_table(str(tmp_path / "nope.fasta"))


def test_load_codon_table_uppercases(tmp_path):
    table_file = tmp_path / "code.fasta"
    table_file.write_text(">atg\nm\n>tgg\nw\n", encoding="utf-8")
    assert attila.load_codon_table(str(table_file)) == {"ATG": "M", "TGG": "W"}


@pytest.mark.parametrize("codon, expected", [
    ("ATG", "M"),
    ("atg", "M"),
    ("NNN", "X"),
    ("AT", ""),
    ("", ""),
])
def test_translate_codon(codon_table, codon, expected):
    assert attila.translate_codon(codon, codon_table) == expected


# --- translate_seq ------------------------------------------------------------

@pytest.mark.parametrize("frame, sense, expected", [
    (1, "+", "MAW"),
    (2, "+", "WPG"),
    (3, "+", "GL"),
    (1, "-", "SRP"),
    (2, "-", "PGH"),
    (3, "-", "QA"),
])
def test_translate_seq_six_frames(codon_table, frame, sense, expected):
    # ATGGCCTGGA, reverse complement TCCAGGCCAT
    assert attila.translate_seq("ATGGCCTGGA", frame, sense, codon_table) == expected


def test_translate_seq_unknown_base_gives_x(codon_table):
    assert attila.translate_seq("ATGNNNTGG", 1, "+", codon_table) == "MXW"


# --- find_cdr3_bounds ---------------------------------------------------------

def _domain(c1_c2_gap, c2_motif_gap, motif="WGQG"):
    """C + gap + C + CDR3 + FR4 motif; Cys1 at 0, Cys2 at gap+1."""
    return "C" + "A" * c1_c2_gap + "C" + "G" * c2_motif_gap + motif


def test_find_cdr3_bounds_vh():
    aa = _domain(55, 7)
    assert attila.find_cdr3_bounds(aa, libtype=0) == (57, 64, 58, 64)
    start, end, _, _ = attila.find_cdr3_bounds(aa, libtype=0)
    assert aa[start:end] == "G" * 7


def test_find_cdr3_bounds_vl_uses_fg_motif():
    vl = _domain(55, 7, motif="FGQG")
    assert attila.find_cdr3_bounds(vl, libtype=1) == (57, 64, 58, 64)
    # The heavy-chain motif is not accepted for light chains and vice versa.
    assert attila.find_cdr3_bounds(_domain(55, 7, motif="WGQG"), libtype=1) is None
    assert attila.find_cdr3_bounds(vl, libtype=0) is None


@pytest.mark.parametrize("gap, found", [(49, False), (50, True), (90, True), (91, False)])
def test_find_cdr3_bounds_cys_distance_limits(gap, found):
    assert (attila.find_cdr3_bounds(_domain(gap, 7), 0) is not None) is found


@pytest.mark.parametrize("gap, found", [(2, False), (3, True), (35, True), (36, False)])
def test_find_cdr3_bounds_cdr3_length_limits(gap, found):
    assert (attila.find_cdr3_bounds(_domain(55, gap), 0) is not None) is found


@pytest.mark.parametrize("aa", [
    "A" * 60 + "WGQG",               # no cysteine
    "C" + "A" * 55 + "C" + "G" * 7,  # no FR4 motif
    "",
])
def test_find_cdr3_bounds_none(aa):
    assert attila.find_cdr3_bounds(aa, 0) is None


# --- translate_read on the synthetic reads ------------------------------------

@pytest.mark.parametrize("variant", sorted(VARIANT_RESIDUE))
@pytest.mark.parametrize("libtype, name, prefix, cdr_prefix, bounds", [
    (0, "vh_final", "VH", "AK", (97, 105)),
    (1, "vl_final", "VL", "", (91, 97)),
])
def test_translate_read_synthetic(codon_table, variant, libtype, name, prefix, cdr_prefix, bounds):
    read_id = f"read_{prefix}_{variant}_1"
    dna = synthetic_read(name, read_id)
    res = attila.translate_read(read_id, dna, codon_table, libtype)
    assert res["frame_str"] == "1+"
    assert res["cdr_seq"] == cdr_prefix + VARIANT_RESIDUE[variant] * 7
    assert (res["cdr_start"], res["cdr_end"]) == bounds
    assert res["aa_seq"] == attila.translate_seq(dna, 1, "+", codon_table)


def test_translate_read_reverse_strand(codon_table):
    dna = attila.reverse_complement(synthetic_read("vh_final", "read_VH_D_1"))
    res = attila.translate_read("r", dna, codon_table, 0)
    assert res["frame_str"] == "1-"
    assert res["cdr_seq"] == "AKQQQQQQQ"


def test_translate_read_shifted_frame(codon_table):
    dna = "AC" + synthetic_read("vh_final", "read_VH_A_1")
    assert attila.translate_read("r", dna, codon_table, 0)["frame_str"] == "3+"


def test_translate_read_without_v_domain(codon_table):
    assert attila.translate_read("junk", "ACGT" * 50, codon_table, 0) is None
