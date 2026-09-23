"""
Run manuscript and supplementary figure scripts, then write paper-named
copies to charts/high_resolution/.

Script names match paper IDs (figure1.py, sup_figure1.py, …). Pair panels
share a primary script (e.g. figure7.py writes Figures 7 and 8).
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

from config import CHARTS_DIR, ensure_charts_dir

ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable

# One entry per generating script (wrappers for the paired panel are omitted).
SCRIPTS = [
    "generate_supplementary_table1.py",
    "figure1.py",
    "figure2.py",
    "figure3.py",
    "figure4a.py",
    "figure4b.py",
    "figure5a.py",
    "figure5b.py",
    "figure6a.py",
    "figure6c.py",
    "figure7.py",
    "figure9.py",
    "figure10.py",
    "figure11.py",
    "figure11m.py",
    "sup_figure1.py",
    "sup_figure2.py",
    "sup_figure3.py",
    "sup_figure5.py",
    "sup_figure6.py",
    "sup_figure8.py",
    "sup_figure10.py",
    "sup_figure12.py",
    "sup_figure14.py",
    "sup_figure16.py",
    "sup_figure18.py",
    "sup_figure19.py",
    "sup_figure20.py",
    "sup_figure21.py",
    "sup_figure22.py",
    "sup_figure23.py",
    "sup_figure24a.py",
    "sup_figure24b.py",
]

POST_SCRIPTS = [
    "export_high_resolution_charts.py",
    "export_figure4b_9b_pixels.py",
]


def _run(script: str, index: int, total: int, env: dict[str, str]) -> bool:
    path = ROOT / script
    if not path.exists():
        print(f"[{index}/{total}] SKIP (not found): {script}")
        return True

    print(f"[{index}/{total}] {script}")
    t0 = time.perf_counter()
    result = subprocess.run([PYTHON, str(path)], cwd=str(ROOT), env=env)
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        print(f"  FAILED (exit {result.returncode}) after {elapsed:.1f}s\n")
        return False
    print(f"  OK ({elapsed:.1f}s)\n")
    return True


def main() -> None:
    ensure_charts_dir()
    print(f"Output folder: {CHARTS_DIR}")
    all_scripts = SCRIPTS + POST_SCRIPTS
    print(f"Running {len(all_scripts)} scripts...\n")

    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    failed: list[str] = []
    total = len(all_scripts)
    for i, script in enumerate(all_scripts, start=1):
        if not _run(script, i, total, env):
            failed.append(script)

    outputs = sorted(p.name for p in CHARTS_DIR.glob("*") if p.is_file())
    hr = CHARTS_DIR / "high_resolution"
    hr_n = len(list(hr.glob("*"))) if hr.exists() else 0
    print(f"Done. {len(outputs)} files in {CHARTS_DIR.name}/; {hr_n} in high_resolution/")
    if failed:
        print(f"Failed scripts ({len(failed)}): {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
