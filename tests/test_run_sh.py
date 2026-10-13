"""Small checks of the run.sh wrapper that do not start the pipeline."""

import re
import shutil
import subprocess

import pytest

from helpers import REPO_ROOT, write

RUN_SH = REPO_ROOT / "run.sh"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")


def _bash_function(name):
    text = RUN_SH.read_text(encoding="utf-8")
    m = re.search(rf"^{name}\(\) \{{\n.*?^\}}\n", text, re.S | re.M)
    assert m, f"{name} not found in run.sh"
    return m.group(0)


def parse_config_value(key, cfg):
    script = _bash_function("parse_config_value") + 'parse_config_value "$1" "$2"\n'
    res = subprocess.run(["bash", "-c", script, "bash", key, str(cfg)], capture_output=True, text=True, check=True)
    return res.stdout.rstrip("\n")


def test_parse_config_value(tmp_path):
    cfg = write(tmp_path / "p_VH.cfg",
                "# projectname: commented out\n"
                "projectname:   my_proj   \n"
                "projectdir:/data/runs\n"
                "igblastdir: /opt/igblast:extra\n")
    assert parse_config_value("projectname", cfg) == "my_proj"
    assert parse_config_value("projectdir", cfg) == "/data/runs"
    assert parse_config_value("igblastdir", cfg) == "/opt/igblast:extra"
    assert parse_config_value("libtype", cfg) == ""


def test_help_exits_zero():
    res = subprocess.run(["bash", str(RUN_SH), "--help"], capture_output=True, text=True, timeout=30)
    assert res.returncode == 0
    assert "Usage: ./run.sh [options]" in res.stdout


def test_unknown_argument_exits_one():
    res = subprocess.run(["bash", str(RUN_SH), "--bogus"], capture_output=True, text=True, timeout=30)
    assert res.returncode == 1
    assert "Unknown argument: --bogus" in res.stdout
