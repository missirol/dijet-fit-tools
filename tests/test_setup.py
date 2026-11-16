"""setup/envrun.sh, setup/activate.sh, setup/shell.sh: a ROOT installation set up in the calling
shell (ROOTSYS, PYTHONPATH, PATH, ...) must not leak into the analysis environment."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MICROMAMBA = os.environ.get("MAMBA_EXE") or shutil.which("micromamba")
pytestmark = [pytest.mark.slow, pytest.mark.skipif(
    not MICROMAMBA or not (Path(sys.prefix) / "conda-meta").is_dir(),
    reason="needs micromamba and a conda environment")]
CHECK = ("import os, ROOT; print('VERSION', ROOT.gROOT.GetVersion()); print('FILE', ROOT.__file__); "
         "print('PYTHONPATH', os.environ.get('PYTHONPATH'))")


@pytest.fixture
def other_root(tmp_path):
    """A fake ROOT installation set up in the calling shell: importing it fails loudly."""
    root = tmp_path / "otherroot"
    (root / "lib").mkdir(parents=True)
    (root / "bin").mkdir()
    (root / "lib" / "ROOT.py").write_text("raise SystemExit('the other ROOT was loaded')\n")
    (root / "bin" / "thisroot.sh").write_text("")
    cfg = root / "bin" / "root-config"
    cfg.write_text("#!/bin/sh\necho 6.00/00\n")
    cfg.chmod(0o755)
    env = dict(os.environ, DIJET_ENV=sys.prefix, MAMBA_EXE=MICROMAMBA, ROOTSYS=str(root),
               PYTHONPATH=f"{root}/lib", LD_LIBRARY_PATH=f"{root}/lib",
               ROOT_INCLUDE_PATH=f"{root}/include", PATH=f"{root}/bin:{os.environ['PATH']}")
    return root, env


def parse(out):
    return dict(line.split(" ", 1) for line in out.splitlines() if " " in line)


def test_envrun_removes_the_other_root(other_root):
    root, env = other_root
    p = subprocess.run([str(REPO / "setup" / "envrun.sh"), "python", "-c", CHECK], env=env,
                       capture_output=True, text=True, cwd="/")
    assert p.returncode == 0, p.stderr
    res = parse(p.stdout)
    assert res["FILE"].startswith(sys.prefix) and res["PYTHONPATH"] == "None"
    assert res["VERSION"] == __import__("ROOT").gROOT.GetVersion()
    assert "ROOTSYS PYTHONPATH LD_LIBRARY_PATH ROOT_INCLUDE_PATH" in p.stderr
    assert f"PATH:{root}/bin" in p.stderr


def test_activate_cleans_and_activates(other_root, tmp_path):
    root, env = other_root
    script = (f"source {REPO}/setup/activate.sh && source {REPO}/setup/activate.sh || exit 7\n"
              f"python -c \"{CHECK}\"\n"
              "echo ROOTCONFIG $(command -v root-config)\n"
              "echo LEFT $(env | grep -c otherroot)\n")
    p = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True, cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    res = parse(p.stdout)
    assert res["FILE"].startswith(sys.prefix) and res["PYTHONPATH"] == "None"
    assert res["ROOTCONFIG"] == f"{sys.prefix}/bin/root-config"
    assert res["LEFT"] == "0"
    assert p.stderr.count("active (ROOT ") == 2


def test_activate_must_be_sourced(other_root):
    _, env = other_root
    p = subprocess.run(["bash", str(REPO / "setup" / "activate.sh")], env=env, capture_output=True,
                       text=True)
    assert p.returncode == 1 and "must be sourced" in p.stderr


def test_shell_reads_bashrc_then_activates(other_root, tmp_path):
    root, env = other_root
    home = tmp_path / "home"
    home.mkdir()
    # the user's startup file sets up the other ROOT again: activate.sh must still win
    (home / ".bashrc").write_text(f"export ROOTSYS={root} PYTHONPATH={root}/lib PATH={root}/bin:$PATH\n"
                                  "alias fromrc='echo FROMRC yes'\n")
    env = dict(env, HOME=str(home))
    p = subprocess.run([str(REPO / "setup" / "shell.sh")], env=env, capture_output=True, text=True,
                       input=f"fromrc\npython -c \"{CHECK}\"\nexit\n", cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    res = parse(p.stdout)
    assert res["FROMRC"] == "yes"
    assert res["FILE"].startswith(sys.prefix) and res["PYTHONPATH"] == "None"


def test_activate_ignores_the_callers_arguments_and_uses_dijet_env(other_root, tmp_path):
    """A sourced script sees the caller's positional parameters: activate.sh takes no arguments
    (regression: they were taken as the environment path). DIJET_ENV selects the environment."""
    _, env = other_root
    script = f"source {REPO}/setup/activate.sh || exit 7\necho PREFIX $CONDA_PREFIX\n"
    p = subprocess.run(["bash", "-c", script, "caller", "not-an-environment"], env=env,
                       capture_output=True, text=True, cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    assert parse(p.stdout)["PREFIX"] == sys.prefix


def test_activate_fails_cleanly_without_an_environment(other_root, tmp_path):
    _, env = other_root
    env = dict(env, DIJET_ENV=str(tmp_path / "nothing"))
    p = subprocess.run(["bash", "-c", f"source {REPO}/setup/activate.sh; echo RC $?"], env=env,
                       capture_output=True, text=True, cwd=tmp_path)
    assert parse(p.stdout)["RC"] != "0" and "no conda environment" in p.stderr


def test_envrun_in_an_active_environment_reports_nothing(other_root, tmp_path):
    """Inside an activated shell the environment's own ROOTSYS and bin/ are not reported as
    'another ROOT installation'."""
    _, env = other_root
    env = {k: v for k, v in env.items() if k not in ("ROOTSYS", "PYTHONPATH", "LD_LIBRARY_PATH",
                                                       "ROOT_INCLUDE_PATH")}
    env["PATH"] = ":".join(d for d in env["PATH"].split(":") if "otherroot" not in d)
    script = (f"source {REPO}/setup/activate.sh 2>/dev/null || exit 7\n"
              f"{REPO}/setup/envrun.sh python -c 'print(1)'\n")
    p = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True, cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    assert "removed from the environment" not in p.stderr, p.stderr


def test_shell_warns_if_the_environment_cannot_be_activated(other_root, tmp_path):
    _, env = other_root
    home = tmp_path / "home"
    home.mkdir()
    env = dict(env, HOME=str(home))
    p = subprocess.run([str(REPO / "setup" / "shell.sh"), str(tmp_path / "nothing")], env=env,
                       capture_output=True, text=True, input="exit\n", cwd=tmp_path)
    assert "NOT active" in p.stderr


def test_clean_env_reports_what_is_not_the_environments_own(tmp_path):
    """Only settings inside the environment (with a '/' boundary) are left out of the report;
    a list that also contains another ROOT, or a sibling directory, is reported."""
    own = tmp_path / "env"
    script = (f". {REPO}/setup/_clean_env.sh\nDIJET_PREFIX={own}\n"
              f"export LD_LIBRARY_PATH={own}/lib:/opt/root/lib ROOTSYS={own}-other PYTHONPATH={own}/lib\n"
              "_dijet_clean_env\necho LEFT ${LD_LIBRARY_PATH:-none} ${ROOTSYS:-none}\n")
    p = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert "LD_LIBRARY_PATH" in p.stderr and "ROOTSYS" in p.stderr and "PYTHONPATH" not in p.stderr
    assert "LEFT none none" in p.stdout
