"""Shared fixtures for the ATTILA test suite.

Nothing here touches the network or runs the real bioinformatics binaries:
`fake_tools` stands in for subprocess.run and `fake_abnum` for the UCL Abnum
web service. Tests marked `tools` opt out of that and use the real programs.
"""

import io
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "programs"))

import attila  # noqa: E402

from helpers import GENCODE, REPO_ROOT, SYNTHETIC_DIR, read_fastq  # noqa: E402


@pytest.fixture(autouse=True)
def _isolate(monkeypatch):
    # attila.main() replaces sys.stderr with an open log file and never restores it.
    monkeypatch.setattr(sys, "stderr", sys.stderr)
    monkeypatch.delenv("ATTILA_STEPS", raising=False)

    def _no_network(*args, **kwargs):
        raise RuntimeError("network access is disabled in tests; use the fake_abnum fixture")

    monkeypatch.setattr(urllib.request, "urlopen", _no_network)


@pytest.fixture
def codon_table():
    return attila.load_codon_table(str(GENCODE))


@pytest.fixture
def make_cfg(tmp_path):
    """Write a .cfg like data/synthetic_data_for_test/example_V?.cfg and return (path, parsed dict)."""

    def _make(libtype=0, paired=False, name=None, **overrides):
        prefix = "vh" if libtype == 0 else "vl"
        values = {
            "projectname": "proj",
            "projectdir": str(tmp_path / "out"),
            "packagedir": str(REPO_ROOT),
            "igblastdir": str(REPO_ROOT / "programs"),
            "libtype": str(libtype),
            "listsize": "5",
            "pairedend": "1" if paired else "0",
            "minlen": "300",
            "minqual": "20",
        }
        if paired:
            values.update({
                "input1r1dir": str(SYNTHETIC_DIR / f"{prefix}_initial.fastq"),
                "input1r2dir": str(SYNTHETIC_DIR / f"{prefix}_initial.fastq"),
                "input2r1dir": str(SYNTHETIC_DIR / f"{prefix}_final.fastq"),
                "input2r2dir": str(SYNTHETIC_DIR / f"{prefix}_final.fastq"),
            })
        else:
            values.update({
                "input1dir": str(SYNTHETIC_DIR / f"{prefix}_initial.fastq"),
                "input2dir": str(SYNTHETIC_DIR / f"{prefix}_final.fastq"),
            })
        values.update({k: str(v) for k, v in overrides.items()})
        lines = ["# Settings for immunoglobulin sequence analysis"]
        lines += [f"{k}: {v}" for k, v in values.items()]
        path = tmp_path / (name or f"proj_{prefix.upper()}.cfg")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path, attila.parse_cfg(str(path))

    return _make


class FakeTools:
    """Stand-in for subprocess.run that records every command and mimics tool outputs."""

    def __init__(self):
        self.calls = []

    def commands(self, tool):
        return [cmd for cmd, _ in self.calls if tool in cmd]

    def __call__(self, cmd, *args, **kwargs):
        cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
        self.calls.append((cmd_str, kwargs))
        tokens = cmd_str.split()
        if tokens and tokens[0] == "prinseq-lite":
            self._prinseq(tokens)
        elif tokens and tokens[0] in ("fastq-join", "fastqjoin"):
            self._fastq_join(tokens)
        elif tokens and tokens[0].endswith("igblastp"):
            self._igblastp(cmd_str, tokens)
        return subprocess.CompletedProcess(cmd, 0, stdout=b"", stderr=b"")

    @staticmethod
    def _opt(tokens, name, default=None):
        return tokens[tokens.index(name) + 1] if name in tokens else default

    def _prinseq(self, tokens):
        reads = read_fastq(self._opt(tokens, "-fastq"))
        min_len = int(self._opt(tokens, "-min_len", 0))
        min_qual = float(self._opt(tokens, "-min_qual_mean", 0))
        good = [
            (rid, seq, qual) for rid, seq, qual in reads
            if len(seq) >= min_len and sum(ord(c) - 33 for c in qual) / max(len(qual), 1) >= min_qual
        ]
        out = Path(self._opt(tokens, "-out_good"))
        out.parent.mkdir(parents=True, exist_ok=True)
        out.with_name(out.name + ".fasta").write_text(
            "".join(f">{rid}\n{seq}\n" for rid, seq, _ in good), encoding="utf-8")
        if self._opt(tokens, "-out_format") == "5":
            out.with_name(out.name + ".fastq").write_text(
                "".join(f"@{rid}\n{seq}\n+\n{qual}\n" for rid, seq, qual in good), encoding="utf-8")

    def _fastq_join(self, tokens):
        # fastq-join appends "join" to the -o argument; the "joined" reads are just R1 here.
        out = Path(self._opt(tokens, "-o") + "join")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(Path(tokens[1]).read_text(encoding="utf-8"), encoding="utf-8")

    def _igblastp(self, cmd_str, tokens):
        dest = Path(cmd_str.split(">")[-1].strip())
        query = Path(self._opt(tokens, "-query"))
        germ = "VH3-23" if "VH_germline" in self._opt(tokens, "-germline_db_V", "") else "V1-13"  # names used in data/germline_human
        lines = []
        for line in query.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                qid = line[1:]
                lines += [
                    f"# Query: {qid}",
                    "FR1\t1\t25", "CDR1\t26\t35", "FR2\t36\t49", "CDR2\t50\t66", "FR3\t67\t98",
                    f"V\t{qid}\t{germ}\t98.2\t100",
                ]
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def fake_tools(monkeypatch):
    fake = FakeTools()
    monkeypatch.setattr(attila.subprocess, "run", fake)
    return fake


class FakeAbnum:
    """Stand-in for urllib.request.urlopen answering like abnum.cgi?plain=1 (one 'H<n> <aa>' per residue)."""

    def __init__(self):
        self.queried = []
        self.fail_for = set()
        self.empty_for = set()

    def __call__(self, req, timeout=None):
        url = req if isinstance(req, str) else req.full_url
        seq = urllib.parse.parse_qs(urllib.parse.urlparse(url).query).get("aaseq", [""])[0]
        self.queried.append(seq)
        if seq in self.fail_for:
            raise OSError("simulated Abnum outage")
        if seq in self.empty_for:
            return io.BytesIO(b"")
        body = "".join(f"H{i} {aa}\n" for i, aa in enumerate(seq, start=1))
        return io.BytesIO(body.encode("utf-8"))


@pytest.fixture
def fake_abnum(monkeypatch):
    fake = FakeAbnum()
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    return fake
