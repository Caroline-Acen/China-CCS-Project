"""
Run all figure scripts sequentially and save outputs to charts/.
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

# Main-text figures first, then supplementary (rough manuscript order)
SCRIPTS = [
    "generate_supplementary_table1.py",
    "figure1.py",
    "figure2.py",
    "figure2b.py",
    "figure3a.py",
    "figure3b.py",
    "figure4a.py",
    "figure4b.py",
    "figure4c.py",
    "figure4d.py",
    "figure5a-b.py",
    "figure6.py",
    "figure6_blue_hydrogen.py",
    "figure7.py",
    "figure1_sup.py",
    "figure2-3_sup.py",
    "figure4-5_sup.py",
    "figure6-7_sup.py",
    "figure8-9_sup.py",
    "figure10-11_sup.py",
    "figure12-13_sup.py",
    "figure14_sup.py",
    "figure14_sup_v2.py",
    "figure15_sup.py",
    "figure16_sup.py",
    "figure17_sup.py",
    "figure18_sup.py",
    "figure19_sup.py",
    "figure20_sup.py",
    "figure21a_sup.py",
    "figure21b_sup.py",
    "figure22_sup.py",
    "figure_ex_en_radial.py",
    "figure_generation_streams_radial.py",
    "figure_weight_sup.py",
    "figure_sink_suitability_distribution.py",
    "run_sink_suitability_map.py",
    "figure_source_sink_cluster_map.py",
    "export_sankey_data.py",
]


def main():
    ensure_charts_dir()
    print(f"Output folder: {CHARTS_DIR}")
    print(f"Running {len(SCRIPTS)} scripts...\n")

    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    failed = []
    for i, script in enumerate(SCRIPTS, start=1):
        path = ROOT / script
        if not path.exists():
            print(f"[{i}/{len(SCRIPTS)}] SKIP (not found): {script}")
            continue

        print(f"[{i}/{len(SCRIPTS)}] {script}")
        t0 = time.perf_counter()
        result = subprocess.run(
            [PYTHON, str(path)],
            cwd=str(ROOT),
            capture_output=False,
            env=env,
        )
        elapsed = time.perf_counter() - t0
        if result.returncode != 0:
            print(f"  FAILED (exit {result.returncode}) after {elapsed:.1f}s\n")
            failed.append(script)
        else:
            print(f"  OK ({elapsed:.1f}s)\n")

    outputs = sorted(CHARTS_DIR.glob("*"))
    print(f"Done. {len(outputs)} files in {CHARTS_DIR.name}/")
    if failed:
        print(f"Failed scripts ({len(failed)}): {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
