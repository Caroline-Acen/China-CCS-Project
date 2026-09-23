"""Export journal-requested pixel sizes for Figure 4b and Figure 9b.

Fonts stay at the original design point sizes (not rescaled with the new canvas).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import Image

from config import CHARTS_DIR, data_path
from figure4b import plot_risk_vs_storage_bubble
from figure9 import (
    compute_hub_metrics,
    generate_blue_hydrogen_gate_maps,
    plot_cost_supply_curve,
)

DPI = 300
HR = CHARTS_DIR / "high_resolution"


def export_figure4b() -> Path:
    # 3150 × 3000 px at 300 dpi; type sized as on the original 14" canvas
    out = HR / "figure4b.png"
    plot_risk_vs_storage_bubble(
        csv_path="./data/Risk_Assessment/final_seismic_risk_factor_base_case.csv",
        output_path=str(out),
        figsize=(3150 / DPI, 3000 / DPI),
        dpi=DPI,
        font_design_width=14.0,
    )
    return out


def export_figure9b() -> Path:
    gasfield_csv = Path(data_path("gasfield_gate_smr_assignments.csv"))
    sink_csv = Path(data_path("sink_gate_smr_assignments.csv"))
    if not (gasfield_csv.exists() and sink_csv.exists()):
        generate_blue_hydrogen_gate_maps()

    gasfield_df = pd.read_csv(gasfield_csv)
    sink_df = pd.read_csv(sink_csv)
    hubs_df = compute_hub_metrics(gasfield_df, sink_df)

    out = HR / "figure9b.png"
    plot_cost_supply_curve(
        hubs_df,
        str(out),
        figsize=(3000 / DPI, 3000 / DPI),
        dpi=DPI,
        font_design_width=10.0,
        font_target_pt=11.0,
        bbox_inches=None,
    )
    return out


def _report(path: Path, expected: tuple[int, int]) -> None:
    with Image.open(path) as im:
        print(f"{path.name}: {im.size[0]}x{im.size[1]} (requested {expected[0]}x{expected[1]})")


if __name__ == "__main__":
    HR.mkdir(parents=True, exist_ok=True)
    p4 = export_figure4b()
    _report(p4, (3150, 3000))
    p9 = export_figure9b()
    _report(p9, (3000, 3000))
