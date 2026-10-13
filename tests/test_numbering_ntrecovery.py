"""Kabat numbering through UCL Abnum (mocked) and recovery of the nucleotide CDS."""

import pytest

import attila
from helpers import read_fasta, synthetic_read, write


ENRICHED = (
    ">read_A|FRAME:1+|FOLD-CHANGE:4.2857|P1:8333333333333.333|P2:35714285714285.72\nEVQLAKGGGWGQG\n"
    ">read_B|FRAME:2-|FOLD-CHANGE:2.0000|P1:1.0|P2:2.0\nQVQLAKSSSWGQG\n"
)


def test_numberab_queries_each_sequence_and_strips_p_values(tmp_path, fake_abnum):
    src = write(tmp_path / "list.fasta", ENRICHED)
    out = tmp_path / "numbered.txt"

    attila.numberab(str(src), str(out))

    assert fake_abnum.queried == ["EVQLAKGGGWGQG", "QVQLAKSSSWGQG"]
    text = out.read_text()
    assert text.startswith(">read_A|FRAME:1+|FOLD-CHANGE:4.2857\nH1 E\nH2 V\n")
    assert ">read_B|FRAME:2-|FOLD-CHANGE:2.0000\n" in text
    assert "P1:" not in text


def test_numberab_skips_failed_and_empty_answers(tmp_path, fake_abnum, capsys):
    src = write(tmp_path / "list.fasta", ENRICHED)
    out = tmp_path / "numbered.txt"
    fake_abnum.fail_for.add("EVQLAKGGGWGQG")
    fake_abnum.empty_for.add("QVQLAKSSSWGQG")

    attila.numberab(str(src), str(out))

    assert out.read_text() == ""
    assert "WARNING: Abnum query failed for read_A" in capsys.readouterr().out


def test_numberab_overwrites_previous_output(tmp_path, fake_abnum):
    src = write(tmp_path / "list.fasta", ENRICHED)
    out = write(tmp_path / "numbered.txt", "stale content\n")
    attila.numberab(str(src), str(out))
    assert "stale" not in out.read_text()


def test_convert_to_fasta_takes_last_residue_letter(tmp_path):
    numbered = write(tmp_path / "numbered.txt",
                     ">id1|FRAME:1+|FOLD-CHANGE:2.0\nH1 E\nH2 V\nH35A -\nH82A Q\nH100 \n"
                     ">id2|FRAME:1+|FOLD-CHANGE:1.5\nL1 D\nL2 I\n")
    out = tmp_path / "numbered.fasta"
    attila.convert_to_fasta(str(numbered), str(out))
    assert read_fasta(out) == [("id1|FRAME:1+|FOLD-CHANGE:2.0", "EVQ"), ("id2|FRAME:1+|FOLD-CHANGE:1.5", "DI")]


def test_numbering_round_trip_reproduces_protein(tmp_path, fake_abnum):
    src = write(tmp_path / "list.fasta", ENRICHED)
    numbered, fasta = tmp_path / "n.txt", tmp_path / "n.fasta"
    attila.numberab(str(src), str(numbered))
    attila.convert_to_fasta(str(numbered), str(fasta))
    assert [seq for _, seq in read_fasta(fasta)] == ["EVQLAKGGGWGQG", "QVQLAKSSSWGQG"]


# --- get_ntsequence -----------------------------------------------------------

def _nt_inputs(tmp_path, codon_table, dna, frame, sense, read_id="read_X"):
    protein = attila.translate_seq(dna, frame, sense, codon_table)
    numbered = write(tmp_path / "numbered.fasta", f">{read_id}|FRAME:{frame}{sense}|FOLD-CHANGE:3.5\n{protein}\n")
    enriched = write(tmp_path / "list.fasta", f">{read_id}|FRAME:{frame}{sense}|FOLD-CHANGE:3.5|P1:1|P2:2\n{protein}\n")
    nt = write(tmp_path / "nt.fasta", f">other|FRAME:1+\nAAAA\n>{read_id}|FRAME:{frame}{sense}\n{dna}\n")
    return numbered, enriched, nt, protein


def test_get_ntsequence_forward_strand(tmp_path, codon_table):
    dna = synthetic_read("vh_final", "read_VH_A_1")
    numbered, enriched, nt, protein = _nt_inputs(tmp_path, codon_table, dna, 1, "+")
    out = tmp_path / "out.fasta"

    attila.get_ntsequence(str(numbered), str(enriched), str(nt), str(out))

    lines = out.read_text().splitlines()
    assert lines[0] == ">read_X|FOLD-CHANGE:3.5"
    assert all(len(l) <= 70 for l in lines[1:])
    assert "".join(lines[1:]) == dna[: len(protein) * 3]


def test_get_ntsequence_reverse_strand(tmp_path, codon_table):
    dna = attila.reverse_complement(synthetic_read("vl_final", "read_VL_E_1"))
    numbered, enriched, nt, protein = _nt_inputs(tmp_path, codon_table, dna, 1, "-")
    out = tmp_path / "out.fasta"

    attila.get_ntsequence(str(numbered), str(enriched), str(nt), str(out))

    cds = "".join(out.read_text().splitlines()[1:])
    # The CDS is returned on the original read's strand.
    assert attila.translate_seq(cds, 1, "-", codon_table) == protein


@pytest.mark.parametrize("protein_has_x, suffix", [(False, ""), (True, " Unknown base (N).")])
def test_get_ntsequence_reports_unrecoverable(tmp_path, protein_has_x, suffix):
    protein = "EVQLXK" if protein_has_x else "EVQLAK"
    numbered = write(tmp_path / "n.fasta", ">id|FRAME:1+|FOLD-CHANGE:2.0\nWWWWWW\n")
    enriched = write(tmp_path / "e.fasta", f">id|FRAME:1+|FOLD-CHANGE:2.0|P1:1|P2:2\n{protein}\n")
    nt = write(tmp_path / "nt.fasta", ">id|FRAME:1+\nGAGGTGCAGCTG\n")
    out = tmp_path / "out.fasta"
    attila.get_ntsequence(str(numbered), str(enriched), str(nt), str(out))
    assert out.read_text() == f">id|FOLD-CHANGE:2.0\nCould not retrieve nucleotide sequence.{suffix}\n"
