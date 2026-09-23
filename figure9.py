"""
Blue hydrogen Dual-Gate SMR map (combined Gasfield-Gate + Sink-Gate).

Hydrogen consumers: MC19 and PM19 only (Blue hydrogen redo.docx).
Distances exclude the 1.3 winding factor by design.

Manuscript figures 9a and 9b: charts/gasfield_sink_combined_gate_map.png
and charts/hydrogen_hub_cost_supply_curve.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from blue_hydrogen_smr import (
    TARGET_CRS,
    gasfield_gate_assignments,
    prepare_blue_hydrogen_layers,
    sink_gate_assignments,
)
from config import (
    BLUE_H2_CO2_TO_H2_RATIO,
    BLUE_H2_H2050_XMAX,
    BLUE_H2_LIFESPAN_YEARS,
    NE_ADMIN0,
    NE_ADMIN1,
    chart_path,
    data_path,
    manuscript_font_bundle,
)

COLOR_CO2_PIPE = "#ff8f00"
COLOR_H2_PIPE = "#c62828"
COLOR_GAS_PIPE = "#6a1b9a"
COLOR_GAS_FIELD = "#d32f2f"
PIPE_LW = 3.8
PIPE_ZORDER = 20
PIPE_ALPHA = 1.0

LOG_PREFIX = "[blue_hydrogen]"


def log(message: str) -> None:
    print(f"{LOG_PREFIX} {message}")


def _load_china_layers():
    admin0 = gpd.read_file(NE_ADMIN0)
    admin1 = gpd.read_file(NE_ADMIN1)
    china = admin0[admin0["ADMIN"] == "China"].to_crs(TARGET_CRS)
    if "admin" in admin1.columns:
        provinces = admin1[admin1["admin"] == "China"]
    elif "ADM0NAME" in admin1.columns:
        provinces = admin1[admin1["ADM0NAME"] == "China"]
    else:
        provinces = gpd.GeoDataFrame(geometry=[], crs=admin1.crs)
    provinces = provinces.to_crs(TARGET_CRS)
    return china, provinces


def _plot_base_map(ax, china, provinces):
    china.boundary.plot(ax=ax, color="black", linewidth=1, zorder=1)
    provinces.boundary.plot(ax=ax, color="gray", linewidth=0.5, alpha=0.7, zorder=1)


def _draw_pipeline(ax, x0, y0, x1, y1, color: str):
    """Draw a high-visibility pipeline with a light halo underneath."""
    ax.plot(
        [x0, x1], [y0, y1],
        color="white",
        linewidth=PIPE_LW + 2.4,
        alpha=0.95,
        zorder=PIPE_ZORDER - 1,
        solid_capstyle="round",
    )
    ax.plot(
        [x0, x1], [y0, y1],
        color=color,
        linewidth=PIPE_LW,
        alpha=PIPE_ALPHA,
        zorder=PIPE_ZORDER,
        solid_capstyle="round",
    )


def _best_sink_gate_spokes(assignments: pd.DataFrame) -> pd.DataFrame:
    """Visible Sink-Gate spokes: best route per gas field, plus longer low-cost routes."""
    if assignments.empty:
        return assignments

    best = assignments.loc[
        assignments.groupby("facility", sort=False)["d_total_weighted"].idxmin()
    ].copy()
    best["leg_km"] = best["d_gas_km"] + best["d_h2_km"]

    tmp = assignments.copy()
    tmp["leg_km"] = tmp["d_gas_km"] + tmp["d_h2_km"]
    low_cost = tmp[tmp["d_total_weighted"] <= tmp["d_total_weighted"].quantile(0.25)]
    long_routes = low_cost.nlargest(40, "leg_km")

    spokes = pd.concat([best, long_routes], ignore_index=True)
    spokes = spokes.drop_duplicates(
        subset=["gas_x", "gas_y", "sink_x", "sink_y", "consumer_x", "consumer_y"]
    )
    return spokes.reset_index(drop=True)


def _add_bottom_legend(fig, handles, legend_rows: int = 3, target_pt: float = 9.0):
    """Place a multi-column legend below the map (no map overlap)."""
    mf = manuscript_font_bundle(fig.get_figwidth(), target_pt=target_pt)
    ncol = max(1, math.ceil(len(handles) / legend_rows))
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.06),
        ncol=ncol,
        fontsize=mf["legend"],
        frameon=True,
        fancybox=False,
        handletextpad=0.4,
        columnspacing=1.2,
        labelspacing=0.5,
        borderaxespad=0.15,
    )
    fig.subplots_adjust(bottom=0.12, top=0.98, left=0.04, right=0.96)


def plot_combined_gate_map(
    gas_gdf,
    sink_gdf,
    consumer_gdf,
    gasfield_assignments: pd.DataFrame,
    sink_assignments: pd.DataFrame,
    output_path: str,
):
    """Single combined map: Gasfield-Gate and Sink-Gate hubs + pipelines."""
    del sink_gdf, consumer_gdf  # used for API compatibility with callers
    log(f"Plotting combined gate map to {output_path}")
    china, provinces = _load_china_layers()
    fig, ax = plt.subplots(figsize=(11, 11))
    _plot_base_map(ax, china, provinces)

    gas_gdf.plot(
        ax=ax, color=COLOR_GAS_FIELD, markersize=48, alpha=0.9,
        edgecolor="black", linewidth=0.45, zorder=8,
    )

    if gasfield_assignments is not None and not gasfield_assignments.empty:
        for row in gasfield_assignments.itertuples():
            _draw_pipeline(ax, row.gas_x, row.gas_y, row.sink_x, row.sink_y, COLOR_CO2_PIPE)
            _draw_pipeline(ax, row.gas_x, row.gas_y, row.consumer_x, row.consumer_y, COLOR_H2_PIPE)

    spokes = (
        _best_sink_gate_spokes(sink_assignments)
        if sink_assignments is not None and not sink_assignments.empty
        else pd.DataFrame()
    )
    if not spokes.empty:
        for row in spokes.itertuples():
            _draw_pipeline(ax, row.gas_x, row.gas_y, row.sink_x, row.sink_y, COLOR_GAS_PIPE)
            _draw_pipeline(ax, row.sink_x, row.sink_y, row.consumer_x, row.consumer_y, COLOR_H2_PIPE)

    if gasfield_assignments is not None and not gasfield_assignments.empty:
        gf_unique = (
            gasfield_assignments[["gas_x", "gas_y"]]
            .drop_duplicates()
            .rename(columns={"gas_x": "x", "gas_y": "y"})
        )
        ax.scatter(
            gf_unique["x"], gf_unique["y"],
            s=180, c="#1e88e5", edgecolors="black", linewidths=0.8, zorder=12,
            label="Gasfield-Gate hub",
        )

    if not spokes.empty:
        sink_unique = (
            spokes[["sink_x", "sink_y"]]
            .drop_duplicates()
            .rename(columns={"sink_x": "x", "sink_y": "y"})
        )
        ax.scatter(
            sink_unique["x"], sink_unique["y"],
            s=140, c="#2e7d32", edgecolors="black", linewidths=0.6, zorder=12,
            label="Sink-Gate hub",
        )

    legend_elements = [
        Patch(facecolor="#1e88e5", edgecolor="black", label="Gasfield-Gate hub"),
        Patch(facecolor="#2e7d32", edgecolor="black", label="Sink-Gate hub"),
        Line2D([0], [0], color=COLOR_CO2_PIPE, linewidth=2.5, label="CO$_2$ pipeline"),
        Line2D([0], [0], color=COLOR_H2_PIPE, linewidth=2.5, label="H$_2$ pipeline"),
        Line2D([0], [0], color=COLOR_GAS_PIPE, linewidth=2.5, label="Natural-gas pipeline"),
    ]
    _add_bottom_legend(fig, legend_elements, legend_rows=2, target_pt=9.0)

    ax.axis("off")
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"Combined Gate map saved to {output_path}")


def compute_hub_metrics(
    gasfield_assignments: pd.DataFrame,
    sink_assignments: pd.DataFrame,
    consumer_gdf: gpd.GeoDataFrame | None = None,
) -> pd.DataFrame:
    """Per-sink hub metrics for the cost-supply curve.

    H_avg,i = S_H2,i / (L * r) with L=30 yr, r=10 (CO2:H2).
    Y_i = min(D_Gasfield-Gate,i, D_Sink-Gate,i) / H_avg,i
    """
    _ = consumer_gdf  # kept for call-site compatibility; H_avg uses storage, not demand
    if sink_assignments is None or sink_assignments.empty:
        log("No sink-gate assignments; cannot compute hub metrics.")
        return pd.DataFrame()

    sk = sink_assignments.copy()
    if "storage_potential" not in sk.columns:
        raise ValueError("sink_assignments must include storage_potential (S_H2)")

    hubs = sk.groupby(["sink_x", "sink_y"], as_index=False).agg(
        D_sink_km=("d_total_weighted", "min"),
        S_H2=("storage_potential", "first"),
        storage_type=("storage_type", "first"),
    )
    hubs["D_sink_km"] = hubs["D_sink_km"] / 1000.0

    if gasfield_assignments is not None and not gasfield_assignments.empty:
        gf = gasfield_assignments.copy()
        gf["D_gasfield_km"] = gf["d_total_weighted"] / 1000.0
        g_by_sink = gf.groupby(["sink_x", "sink_y"], as_index=False)["D_gasfield_km"].min()
        hubs = hubs.merge(g_by_sink, on=["sink_x", "sink_y"], how="left")
    else:
        hubs["D_gasfield_km"] = np.nan

    denom = BLUE_H2_LIFESPAN_YEARS * BLUE_H2_CO2_TO_H2_RATIO
    hubs["H_avg"] = hubs["S_H2"] / denom
    hubs["D_min_km"] = hubs[["D_gasfield_km", "D_sink_km"]].min(axis=1, skipna=True)
    hubs["Y"] = hubs["D_min_km"] / hubs["H_avg"].replace({0: np.nan})
    hubs = hubs.rename(columns={"sink_x": "x", "sink_y": "y"})
    hubs["hub_type"] = np.where(
        hubs["D_gasfield_km"].notna() & (hubs["D_gasfield_km"] <= hubs["D_sink_km"]),
        "gasfield",
        "sink",
    )

    log(f"Computed hubs: {len(hubs)} rows")
    return hubs


def plot_cost_supply_curve(
    hubs_df: pd.DataFrame,
    output_path: str,
    *,
    figsize=(6.5, 4.2),
    dpi=300,
    font_design_width=None,
    font_target_pt=9.0,
    bbox_inches="tight",
):
    """Stepwise (MAC-style) cost-supply curve: Y vs cumulative H_avg up to H_2050 max."""
    if hubs_df is None or hubs_df.empty:
        log("No hubs available to plot cost-supply curve.")
        return

    df = hubs_df.dropna(subset=["Y", "H_avg"]).copy()
    df = df[df["H_avg"] > 0]
    if df.empty:
        log("No hubs with positive H_avg to plot.")
        return

    df_sorted = df.sort_values("Y").reset_index(drop=True)
    xmax = float(BLUE_H2_H2050_XMAX)

    # Keep cheapest hubs until cumulative capacity reaches the Table-4 H_2050 cap
    cum = 0.0
    widths: list[float] = []
    heights: list[float] = []
    lefts: list[float] = []
    for _, row in df_sorted.iterrows():
        if cum >= xmax:
            break
        width = float(row["H_avg"])
        if cum + width > xmax:
            width = xmax - cum
        lefts.append(cum)
        widths.append(width)
        heights.append(float(row["Y"]))
        cum += width

    if not widths:
        log("No hubs fall within the H_2050 x-axis range.")
        return

    fig, ax = plt.subplots(figsize=figsize)
    mf = manuscript_font_bundle(
        font_design_width if font_design_width is not None else fig.get_figwidth(),
        target_pt=font_target_pt,
    )
    ax.bar(
        lefts,
        heights,
        width=widths,
        align="edge",
        color="#1e88e5",
        edgecolor="#0d47a1",
        linewidth=0.6,
    )
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, max(heights) * 1.08)
    ax.set_xlabel(
        r"Cumulative Blue Hydrogen Production Capacity by 2050 (Mt H$_2$/year)",
        fontsize=mf["label"],
    )
    ax.set_ylabel("Levelized Network Friction Index (Y)", fontsize=mf["label"])
    ax.tick_params(labelsize=mf["tick"])
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(
        output_path, dpi=dpi, bbox_inches=bbox_inches, facecolor="white", pad_inches=0
    )
    plt.close(fig)
    log(f"Cost-supply curve saved to {output_path} (x=0–{xmax:g}, {len(widths)} steps)")


def generate_blue_hydrogen_gate_maps(
    *, reuse_assignments: bool | None = None, recompute_assignments: bool = False
):
    log("Loading blue hydrogen layers (MC19 + PM19)...")
    gas_gdf, sink_gdf, consumer_gdf = prepare_blue_hydrogen_layers()
    log(
        f"Gas fields: {len(gas_gdf):,}; storage sites: {len(sink_gdf):,}; "
        f"H2 consumers: {len(consumer_gdf):,}"
    )

    gasfield_csv = Path(data_path("gasfield_gate_smr_assignments.csv"))
    sink_csv = Path(data_path("sink_gate_smr_assignments.csv"))
    have_saved = gasfield_csv.exists() and sink_csv.exists()
    reuse = (not recompute_assignments) and (
        True if reuse_assignments is None else reuse_assignments
    ) and have_saved

    if reuse:
        log("Reusing saved assignment CSVs...")
        gasfield_df = pd.read_csv(gasfield_csv)
        sink_df = pd.read_csv(sink_csv)
    else:
        log("Computing Gasfield-Gate assignments...")
        gasfield_df = gasfield_gate_assignments(gas_gdf, sink_gdf, consumer_gdf)
        gasfield_df.to_csv(gasfield_csv, index=False)
        log(f"{len(gasfield_df):,} viable gas-field gates -> {gasfield_csv}")

        log("Computing Sink-Gate assignments...")
        sink_df = sink_gate_assignments(gas_gdf, sink_gdf, consumer_gdf)
        sink_df.to_csv(sink_csv, index=False)
        log(f"{len(sink_df):,} viable sink gates -> {sink_csv}")

    log(f"Gasfield-Gate pairs: {len(gasfield_df):,}; Sink-Gate pairs: {len(sink_df):,}")
    plot_combined_gate_map(
        gas_gdf,
        sink_gdf,
        consumer_gdf,
        gasfield_df,
        sink_df,
        chart_path("gasfield_sink_combined_gate_map.png"),
    )

    hubs_df = compute_hub_metrics(gasfield_df, sink_df, consumer_gdf)
    plot_cost_supply_curve(hubs_df, chart_path("hydrogen_hub_cost_supply_curve.png"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Blue hydrogen Dual-Gate combined SMR map")
    parser.add_argument(
        "--reuse-assignments",
        action="store_true",
        help="Reuse saved assignment CSVs (default when those files exist)",
    )
    parser.add_argument(
        "--recompute-assignments",
        action="store_true",
        help="Recompute Dual-Gate assignments even if CSVs already exist",
    )
    args = parser.parse_args()
    generate_blue_hydrogen_gate_maps(
        reuse_assignments=True if args.reuse_assignments else None,
        recompute_assignments=args.recompute_assignments,
    )
