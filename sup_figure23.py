"""Supplementary Figure 23: Weight-sensitivity SRF maps."""

Scenarios from Risk Assessment Redo.docx (last table):
  (a) Tectonic Dominance       w = (0.30, 0.50, 0.20)
  (b) Geotechnical Dominance   w = (0.40, 0.20, 0.40)
  (c) Absolute Shaking         w = (0.60, 0.30, 0.10)

Reweights SPGA / SFP / SSA from the base-case CSV (no full recompute).
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from sklearn.neighbors import NearestNeighbors

from config import NE_ADMIN0, NE_ADMIN1, RISK_DIR, chart_path, manuscript_font_bundle
from seismic_risk import (
    SSI_COLORBAR_LABEL,
    SSI_RISK_LEGEND,
    WEIGHT_SCENARIOS,
    reweight_risk_table,
    srf_listed_colormap,
    style_srf_colorbar,
)

EXPORT_COLS = [
    "lat", "lon", "storage_type", "storage_potential", "final_risk_score",
    "pga_level", "fault_lambda", "amplification_method",
    "pga_score", "fault_score", "site_vulnerability",
    "analysis_type", "risk_category", "w1", "w2", "w3",
]


def write_weight_csvs(base_csv: str | None = None) -> list[tuple[str, str]]:
    """Return list of (csv_path, panel_title) for the three weight scenarios."""
    base_path = RISK_DIR / (base_csv or "final_seismic_risk_factor_base_case.csv")
    base = pd.read_csv(base_path)
    required = {"pga_score", "fault_score", "site_vulnerability"}
    missing = required - set(base.columns)
    if missing:
        raise ValueError(f"Base CSV missing component columns: {missing}")

    out_meta: list[tuple[str, str]] = []
    frames: list[pd.DataFrame] = []
    for filename, title, (w1, w2, w3) in WEIGHT_SCENARIOS:
        df = reweight_risk_table(base, w1, w2, w3)
        path = RISK_DIR / filename
        df[EXPORT_COLS].to_csv(path, index=False)
        print(f"[OK] {path.name}  w=({w1:.2f}, {w2:.2f}, {w3:.2f})  n={len(df):,}")
        out_meta.append((str(path), title))
        frames.append(df.assign(panel="abc"[len(frames)], scenario=title))

    # Combined wide CSV for the figure (one row per site, SRF for a/b/c)
    a, b, c = frames
    key = [
        "lat", "lon", "storage_type", "storage_potential",
        "pga_score", "fault_score", "site_vulnerability",
    ]
    wide = a[key].copy()
    wide["srf_a_tectonic"] = a["final_risk_score"]
    wide["risk_category_a"] = a["risk_category"]
    wide["srf_b_geotechnical"] = b["final_risk_score"]
    wide["risk_category_b"] = b["risk_category"]
    wide["srf_c_shaking"] = c["final_risk_score"]
    wide["risk_category_c"] = c["risk_category"]
    wide_path = RISK_DIR / "sensitivity_weight_abc.csv"
    wide.to_csv(wide_path, index=False)
    print(f"[OK] combined wide CSV -> {wide_path}")

    long = pd.concat(frames, ignore_index=True)[
        [
            "panel", "scenario", "w1", "w2", "w3",
            "lat", "lon", "storage_type", "storage_potential",
            "pga_score", "fault_score", "site_vulnerability",
            "final_risk_score", "risk_category",
        ]
    ]
    long_path = RISK_DIR / "sensitivity_weight_abc_long.csv"
    long.to_csv(long_path, index=False)
    print(f"[OK] combined long CSV -> {long_path}")

    return out_meta


def plot_weight_sensitivity(
    panels: list[tuple[str, str]],
    admin0_path: str,
    admin1_path: str,
    output_path: str,
    target_crs: str = "ESRI:102012",
    link_dist_m: float = 60_000,
    k_neighbors: int = 4,
    line_width: float = 0.18,
    line_alpha: float = 1.0,
    point_size: float = 1.2,
    point_alpha: float = 1.0,
    dpi: int = 600,
) -> None:
    n_panels = len(panels)
    ncols, nrows = n_panels, 1

    admin0 = gpd.read_file(admin0_path).to_crs(target_crs)
    admin1 = gpd.read_file(admin1_path).to_crs(target_crs)

    china0 = None
    for col in ["ADMIN", "NAME", "SOVEREIGNT", "ADM0NAME", "NAME_EN", "COUNTRY"]:
        if col in admin0.columns:
            china0 = admin0[admin0[col].astype(str).str.strip().str.lower().eq("china")]
            if not china0.empty:
                break
    if china0 is None or china0.empty:
        raise ValueError("China not found in admin0 shapefile.")

    provinces = None
    for col in ["admin", "ADM0NAME", "ADMIN", "SOVEREIGNT", "COUNTRY"]:
        if col in admin1.columns:
            provinces = admin1[admin1[col].astype(str).str.strip().str.lower().eq("china")]
            if provinces is not None and not provinces.empty:
                break
    if provinces is None or provinces.empty:
        provinces = admin1.copy()

    try:
        provinces = gpd.overlay(provinces, china0[["geometry"]], how="intersection")
    except Exception:
        provinces = provinces.copy()
        provinces["geometry"] = provinces.geometry.buffer(0)
        china_fixed = china0.copy()
        china_fixed["geometry"] = china_fixed.geometry.buffer(0)
        provinces = gpd.overlay(provinces, china_fixed[["geometry"]], how="intersection")

    minx, miny, maxx, maxy = china0.total_bounds
    pad_x = (maxx - minx) * 0.02
    pad_y = (maxy - miny) * 0.02
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y)

    gdfs: list[gpd.GeoDataFrame] = []
    global_max = 0.0
    for csv_path, _title in panels:
        df = pd.read_csv(csv_path)
        df = df[["lat", "lon", "final_risk_score"]].dropna()
        df["final_risk_score"] = pd.to_numeric(df["final_risk_score"], errors="coerce")
        df = df.dropna(subset=["final_risk_score"]).drop_duplicates()
        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(df["lon"], df["lat"]),
            crs="EPSG:4326",
        ).to_crs(target_crs)
        gdfs.append(gdf)
        if len(gdf):
            global_max = max(global_max, float(gdf["final_risk_score"].max()))

    vmin = 0
    vmax = int(np.ceil(global_max))
    cmap, norm, boundaries = srf_listed_colormap(vmin, vmax)

    if vmax <= 12:
        step = 1
    elif vmax <= 25:
        step = 2
    elif vmax <= 50:
        step = 5
    else:
        step = 10
    ticks = list(np.arange(vmin, vmax + 1, step))
    if ticks[-1] != vmax:
        ticks.append(vmax)

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(11 * ncols, 7 * nrows), constrained_layout=False,
    )
    mf = manuscript_font_bundle(fig.get_figwidth())
    axes = np.array(axes).reshape(-1)
    # Match amplification sensitivity layout (same per-map size + margins)
    fig.subplots_adjust(left=0.02, right=0.80, top=0.98, bottom=0.22, wspace=0.02, hspace=0.10)

    letters = "abcdefghijklmnopqrstuvwxyz"
    for i, ax in enumerate(axes):
        if i >= n_panels:
            ax.axis("off")
            continue

        gdf = gdfs[i]
        china0.boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
        provinces.boundary.plot(ax=ax, color="0.75", linewidth=0.6, zorder=2)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.axis("off")

        lab = f"({letters[i]})" if i < len(letters) else f"({i+1})"
        ax.text(
            0.5, 0.995, lab, transform=ax.transAxes,
            ha="center", va="top", fontsize=mf["title"], fontweight="bold", zorder=10,
        )

        xy = np.column_stack([gdf.geometry.x.values, gdf.geometry.y.values])
        scores = gdf["final_risk_score"].values
        segments, seg_vals = [], []
        n = len(xy)
        if n >= 2:
            k = min(int(k_neighbors), n)
            nn = NearestNeighbors(n_neighbors=k, algorithm="auto")
            nn.fit(xy)
            dists, idxs = nn.kneighbors(xy)
            for p in range(n):
                for jpos in range(1, k):
                    q = idxs[p, jpos]
                    if dists[p, jpos] <= link_dist_m:
                        segments.append([xy[p], xy[q]])
                        seg_vals.append(scores[p])

        if segments:
            lc = LineCollection(
                segments, cmap=cmap, norm=norm, linewidths=line_width,
                alpha=line_alpha, zorder=3, rasterized=True,
            )
            lc.set_array(np.asarray(seg_vals))
            ax.add_collection(lc)

        ax.scatter(
            xy[:, 0], xy[:, 1], c=scores, cmap=cmap, norm=norm,
            s=point_size, linewidths=0, alpha=point_alpha, zorder=4, rasterized=True,
        )

    # Stretch colorbar to the full height of the map axes (like amplification figure)
    active = [axes[i] for i in range(n_panels)]
    y0 = min(ax.get_position().y0 for ax in active)
    y1 = max(ax.get_position().y1 for ax in active)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.83, y0, 0.015, y1 - y0])
    cbar = plt.colorbar(
        sm, cax=cax, orientation="vertical",
        boundaries=boundaries, ticks=ticks, spacing="proportional", drawedges=True,
    )
    style_srf_colorbar(
        fig,
        cbar,
        SSI_COLORBAR_LABEL,
        mf["colorbar_label"],
        round(mf["colorbar_tick"] * 0.95, 1),
        label_dx=0.04,
    )

    fig.text(0.5, 0.055, SSI_RISK_LEGEND, ha="center", va="center", fontsize=mf["legend"])

    plt.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print(f"[OK] Saved {output_path}")


if __name__ == "__main__":
    panels = write_weight_csvs()
    plot_weight_sensitivity(
        panels=panels,
        admin0_path=str(NE_ADMIN0),
        admin1_path=str(NE_ADMIN1),
        output_path=chart_path("sensitivity_weight_abc.png"),
    )
