"""
Recompute seismic risk CSVs (updated directional SFP) and optionally replot maps.

Example:
    python compute_seismic_risk.py --scenario base
    python compute_seismic_risk.py --all
    python compute_seismic_risk.py --all --plot
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from config import RISK_DIR, ensure_charts_dir
from seismic_risk import DEFAULT_SCENARIOS, RiskScenario, compute_risk_table

EXPORT_COLS = [
    "lat", "lon", "storage_type", "storage_potential", "final_risk_score",
    "pga_level", "fault_lambda", "amplification_method",
    "pga_score", "fault_score", "site_vulnerability",
    "analysis_type", "risk_category",
]


def _write_scenario(filename: str, scenario: RiskScenario) -> Path:
    out = RISK_DIR / filename
    t0 = time.perf_counter()
    df = compute_risk_table(scenario)
    export = df.rename(columns={"lat": "lat", "lon": "lon"})
    export[EXPORT_COLS].to_csv(out, index=False)
    elapsed = time.perf_counter() - t0
    print(f"[OK] {out.name} ({len(export):,} rows, {elapsed:.1f}s)")
    return out


def _plot_figures() -> None:
    scripts = [
        "figure3a.py",
        "figure3b.py",
        "figure17_sup.py",
        "figure18_sup.py",
        "figure19_sup.py",
    ]
    root = Path(__file__).resolve().parent
    ensure_charts_dir()
    for name in scripts:
        path = root / name
        if not path.exists():
            continue
        print(f"[plot] {name}")
        subprocess.run([sys.executable, str(path)], check=True, cwd=str(root))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Recompute seismic risk assessment CSVs.")
    p.add_argument("--all", action="store_true", help="Write all sensitivity + base-case CSVs.")
    p.add_argument(
        "--scenario",
        choices=("base", "pga", "lambda", "amplification"),
        help="Write one scenario group.",
    )
    p.add_argument("--plot", action="store_true", help="Run risk map figure scripts after compute.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.all and not args.scenario:
        args.scenario = "base"

    selected: list[tuple[str, RiskScenario]] = []
    if args.all:
        selected = DEFAULT_SCENARIOS
    elif args.scenario == "base":
        selected = [DEFAULT_SCENARIOS[0]]
    elif args.scenario == "pga":
        selected = [s for s in DEFAULT_SCENARIOS if "pga" in s[0]]
    elif args.scenario == "lambda":
        selected = [s for s in DEFAULT_SCENARIOS if "lambda" in s[0]]
    elif args.scenario == "amplification":
        selected = [s for s in DEFAULT_SCENARIOS if "amplification" in s[0]]

    for filename, scenario in selected:
        _write_scenario(filename, scenario)

    if args.plot:
        _plot_figures()


if __name__ == "__main__":
    main()
