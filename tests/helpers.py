"""Small file helpers shared by the test modules."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PROGRAMS_DIR = REPO_ROOT / "programs"
SYNTHETIC_DIR = REPO_ROOT / "data" / "synthetic_data_for_test"
GENCODE = REPO_ROOT / "data" / "genetic_code.fasta"

# Synthetic dataset design (see data/synthetic_data_for_test/generate_synthetic_data.py):
# each variant is a CDR3 of 7 identical residues; reads per round are R0 -> RN.
VARIANT_RESIDUE = {"A": "G", "B": "S", "C": "A", "D": "Q", "E": "Y"}
R0_COUNTS = {"A": 1, "B": 5, "C": 3, "D": 1, "E": 2}   # 12 reads
RN_COUNTS = {"A": 10, "B": 1, "C": 3, "D": 8, "E": 6}  # 28 reads


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def read_fasta(path):
    """Return a list of (header_without_gt, sequence) tuples."""
    records = []
    header, seq = None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        elif header is not None:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def read_fastq(path):
    """Return a list of (id, sequence, quality) tuples."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines) - 3, 4)]


def synthetic_read(name, read_id):
    """DNA of one read from a synthetic FASTQ, e.g. synthetic_read('vh_final', 'read_VH_A_1')."""
    for rid, seq, _ in read_fastq(SYNTHETIC_DIR / f"{name}.fastq"):
        if rid == read_id:
            return seq
    raise KeyError(read_id)


def lib_dirs(cfg):
    """Create and return the InitialRound/FinalRound/SelectedSequences/Report dirs for a cfg dict."""
    lib = "VH" if int(cfg["libtype"]) == 0 else "VL"
    base = Path(cfg["projectdir"]) / cfg["projectname"]
    dirs = {
        "lib": base / lib,
        "init": base / lib / "InitialRound",
        "final": base / lib / "FinalRound",
        "selected": base / lib / "SelectedSequences",
        "report": base / "Report",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs
