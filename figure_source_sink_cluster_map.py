"""
China source–sink cluster map (complement to sink_suitability_pareto.png).

Full-China panels by year: provincial basemap, top-15% SS_ij flows (orange lines),
industry context, DSA/EOR sinks. No figure suptitle.

Example:
    python figure_source_sink_cluster_map.py --categories MC
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from shapely.geometry import LineString, Point, box as shapely_box
from sklearn.neighbors import BallTree

from config import COST_YEARS, NE_ADMIN0, NE_ADMIN1, SS_SCORES_DIR, chart_path, ensure_charts_dir, manuscript_font_bundle
from figure_sink_suitability_distribution import (
    INDUSTRY_CATEGORIES,
    load_sinks_for_year,
    _max_ss_for_emitters,
)

# All cost years in one national figure (2030 uses the same pipeline as other years)
DEFAULT_MAP_YEARS = list(COST_YEARS)

CHINA_ALBERS = "EPSG:4479"

HUB_LINE_COLORS = ("#9467bd", "#2ca02c", "#d62728", "#1f77b4")
COLOR_HUB_MARKER = "#d62728"

# Overview map symbology
COLOR_INDUSTRY_BG = "#9e9e9e"
COLOR_INDUSTRY_MATCH = "#08519c"
COLOR_FLOW = "#ff6600"
COLOR_FLOW_HALO = "#ffffff"
COLOR_DSA = "#1b7837"
COLOR_EOR = "#762a83"

# Regional single-panel symbology
COLOR_SOURCE_REGIONAL = "#2ca02c"
COLOR_FLOW_REGIONAL = "#c0392b"

MATCH_REGIONS: dict[str, tuple[float, float, float, float]] = {
    "East": (116.5, 121.0, 29.5, 34.5),
    "North-China": (115.5, 120.0, 36.0, 40.5),
    "Northeast": (121.5, 127.5, 41.0, 46.0),
    "West": (75.5, 82.5, 37.5, 42.0),
}


def _resolve_categories(names: list[str] | None) -> list[str]:
    if not names:
        return list(INDUSTRY_CATEGORIES.values())
    out = []
    for n in names:
        n = n.strip().upper()
        if n in INDUSTRY_CATEGORIES:
            out.append(INDUSTRY_CATEGORIES[n])
        elif n in INDUSTRY_CATEGORIES.values():
            out.append(n)
        else:
            raise ValueError(f"Unknown category {n!r}")
    return out


def _enrich_sink_coords(df: pd.DataFrame, year: int) -> pd.DataFrame:
    if df.empty or {"Best_Sink_X", "Best_Sink_Y"}.issubset(df.columns):
        return df
    sinks = load_sinks_for_year(year, storage_mode="all")
    sink_x = sinks["X"].to_numpy(dtype=np.float64)
    sink_y = sinks["Y"].to_numpy(dtype=np.float64)
    sink_storage_type = sinks["StorageType"].to_numpy()
    r_t = sinks["R_t"].to_numpy()
    s_t = sinks["S_t"].to_numpy()
    tree = BallTree(
        np.column_stack([np.radians(sink_y), np.radians(sink_x)]),
        metric="haversine",
    )
    lat = df["Industry_Y"].to_numpy(dtype=np.float64)
    lon = df["Industry_X"].to_numpy(dtype=np.float64)
    _, best_dist, best_si = _max_ss_for_emitters(lat, lon, tree, r_t, s_t)
    out = df.copy()
    out["Best_Sink_X"] = sink_x[best_si]
    out["Best_Sink_Y"] = sink_y[best_si]
    out["Sink_StorageType"] = sink_storage_type[best_si]
    if "Best_Distance_km" not in out.columns:
        out["Best_Distance_km"] = best_dist
    return out


def _thin_matches_from_files(
    year: int,
    categories: list[str],
    ss_quantile: float,
    max_lines: int,
    *,
    min_link_km: float = 0.0,
) -> pd.DataFrame:
    """
    Keep rows in the top (1 - ss_quantile) fraction by Max_SS (same rule as Pareto).

    When min_link_km > 0, prefer separable source→sink pairs for map lines. Taking only
    nlargest(Max_SS) would keep co-located pairs (SS ∝ 1/d²) and hide all spokes.
    """
    sink_cols = ["Best_Sink_X", "Best_Sink_Y", "Sink_StorageType"]
    frames: list[pd.DataFrame] = []

    for short in categories:
        p = SS_SCORES_DIR / f"ss_{short}_{year}.csv"
        if not p.exists():
            continue
        header = pd.read_csv(p, nrows=0).columns.tolist()
        has_sink = all(c in header for c in sink_cols)
        usecols = [
            "Industry_X", "Industry_Y", "Max_SS", "Category", "Best_Distance_km",
            *(sink_cols if has_sink else []),
        ]
        usecols = [c for c in usecols if c in header]
        chunks_top: list[pd.DataFrame] = []
        for chunk in pd.read_csv(p, usecols=usecols, chunksize=200_000):
            chunk = chunk[chunk["Max_SS"] > 0]
            if not chunk.empty:
                chunks_top.append(chunk)
        if chunks_top:
            frames.append(pd.concat(chunks_top, ignore_index=True))

    if not frames:
        raise FileNotFoundError(f"No score data for year={year}.")

    pool = pd.concat(frames, ignore_index=True)
    q = pool["Max_SS"].quantile(ss_quantile)
    high = pool[pool["Max_SS"] >= q]
    if high.empty:
        return high

    if min_link_km > 0 and "Best_Distance_km" in high.columns:
        linked = high[high["Best_Distance_km"] >= min_link_km]
        if len(linked) >= max_lines:
            return linked.nlargest(max_lines, "Max_SS")
        chosen = linked.nlargest(len(linked), "Max_SS")
        need = max_lines - len(chosen)
        if need > 0:
            short = high[high["Best_Distance_km"] < min_link_km]
            if not short.empty:
                extra = short.nlargest(min(need, len(short)), "Max_SS")
                chosen = pd.concat([chosen, extra], ignore_index=True)
        return chosen

    if len(high) > max_lines:
        return high.sample(n=max_lines, random_state=year)
    return high


def _sample_industries(year: int, categories: list[str], max_points: int, *, seed: int) -> pd.DataFrame:
    usecols = ["Industry_X", "Industry_Y", "Max_SS", "Category"]
    frames: list[pd.DataFrame] = []
    n_paths = sum(1 for s in categories if (SS_SCORES_DIR / f"ss_{s}_{year}.csv").exists())
    per_file = max(2000, max_points // max(1, n_paths))

    for short in categories:
        p = SS_SCORES_DIR / f"ss_{short}_{year}.csv"
        if not p.exists():
            continue
        taken = 0
        for chunk_idx, chunk in enumerate(pd.read_csv(p, usecols=usecols, chunksize=400_000)):
            chunk = chunk[chunk["Max_SS"] > 0]
            if chunk.empty:
                continue
            n_draw = min(len(chunk), per_file // 4 + 800)
            frames.append(chunk.sample(n=n_draw, random_state=seed + chunk_idx + year))
            taken += n_draw
            if taken >= per_file:
                break

    if not frames:
        return pd.DataFrame(columns=["Industry_X", "Industry_Y"])
    out = pd.concat(frames, ignore_index=True)
    if len(out) > max_points:
        out = out.sample(n=max_points, random_state=seed)
    return out[["Industry_X", "Industry_Y"]]


def _in_box(df: pd.DataFrame, box: tuple[float, float, float, float]) -> pd.Series:
    x0, x1, y0, y1 = box
    ind = df["Industry_X"].between(x0, x1) & df["Industry_Y"].between(y0, y1)
    if {"Best_Sink_X", "Best_Sink_Y"}.issubset(df.columns):
        sink = df["Best_Sink_X"].between(x0, x1) & df["Best_Sink_Y"].between(y0, y1)
        return ind | sink
    return ind


def _region_for_row(lon: float, lat: float) -> str | None:
    for name, (x0, x1, y0, y1) in MATCH_REGIONS.items():
        if x0 <= lon <= x1 and y0 <= lat <= y1:
            return name
    return None


def _load_best_region_pool(
    year: int,
    categories: list[str],
    min_link_km: float,
    pool_size: int,
) -> tuple[pd.DataFrame, str]:
    """One CSV pass per year; assign rows to macro-regions; pick the richest."""
    sink_cols = ["Best_Sink_X", "Best_Sink_Y", "Sink_StorageType"]
    buckets: dict[str, list[pd.DataFrame]] = {k: [] for k in MATCH_REGIONS}

    for short in categories:
        p = SS_SCORES_DIR / f"ss_{short}_{year}.csv"
        if not p.exists():
            continue
        header = pd.read_csv(p, nrows=0).columns.tolist()
        has_sink = all(c in header for c in sink_cols)
        usecols = [
            "Industry_X", "Industry_Y", "Max_SS", "Category", "Best_Distance_km",
            *(sink_cols if has_sink else []),
        ]
        usecols = [c for c in usecols if c in header]
        for chunk in pd.read_csv(p, usecols=usecols, chunksize=400_000):
            chunk = chunk[chunk["Max_SS"] > 0]
            if "Best_Distance_km" in chunk.columns:
                chunk = chunk[chunk["Best_Distance_km"] >= min_link_km]
            if chunk.empty:
                continue
            for name, box in MATCH_REGIONS.items():
                sub = chunk[_in_box(chunk, box)]
                if not sub.empty:
                    buckets[name].append(sub)

    best_name = "China"
    best_df = pd.DataFrame()
    for name, parts in buckets.items():
        if not parts:
            continue
        merged = pd.concat(parts, ignore_index=True)
        merged = _enrich_sink_coords(merged, year)
        merged = merged[merged["Best_Distance_km"] >= min_link_km]
        merged = merged.nlargest(min(pool_size, len(merged)), "Max_SS")
        if len(merged) > len(best_df):
            best_name, best_df = name, merged

    if best_df.empty:
        raise FileNotFoundError(f"No links ≥ {min_link_km} km for year={year}.")
    return best_df, best_name


def _load_spoke_pool_in_box(
    year: int,
    categories: list[str],
    box: tuple[float, float, float, float],
    min_link_km: float,
    pool_size: int,
) -> pd.DataFrame:
    """Long-link candidates within a macro-region window."""
    sink_cols = ["Best_Sink_X", "Best_Sink_Y", "Sink_StorageType"]
    parts: list[pd.DataFrame] = []

    for short in categories:
        p = SS_SCORES_DIR / f"ss_{short}_{year}.csv"
        if not p.exists():
            continue
        header = pd.read_csv(p, nrows=0).columns.tolist()
        has_sink = all(c in header for c in sink_cols)
        usecols = [
            "Industry_X", "Industry_Y", "Max_SS", "Category", "Best_Distance_km",
            *(sink_cols if has_sink else []),
        ]
        usecols = [c for c in usecols if c in header]
        for chunk in pd.read_csv(p, usecols=usecols, chunksize=300_000):
            chunk = chunk[chunk["Max_SS"] > 0]
            if "Best_Distance_km" in chunk.columns:
                chunk = chunk[chunk["Best_Distance_km"] >= min_link_km]
            chunk = chunk[_in_box(chunk, box)]
            if not chunk.empty:
                parts.append(chunk)

    if not parts:
        return pd.DataFrame()
    merged = pd.concat(parts, ignore_index=True)
    merged = _enrich_sink_coords(merged, year)
    merged = merged[merged["Best_Distance_km"] >= min_link_km]
    return merged.nlargest(min(pool_size, len(merged)), "Max_SS")


def _load_spoke_pool(
    year: int,
    categories: list[str],
    min_link_km: float,
    pool_size: int,
) -> pd.DataFrame:
    df, _ = _load_best_region_pool(year, categories, min_link_km, pool_size)
    return df


def _select_hub_spokes(
    df: pd.DataFrame,
    *,
    ss_quantile: float,
    min_link_km: float,
    max_hubs: int,
    max_per_hub: int,
) -> pd.DataFrame:
    """
    Top SS_ij pairs with visible separation, grouped by storage hub (sink).
    Mirrors Pareto top-15% filter; keeps only links long enough to draw as spokes.
    """
    if df.empty:
        return df
    sub = df.copy()
    sub = sub[
        (sub["Industry_X"] != sub["Best_Sink_X"]) | (sub["Industry_Y"] != sub["Best_Sink_Y"])
    ]
    if "Best_Distance_km" in sub.columns:
        sub = sub[sub["Best_Distance_km"] >= min_link_km]
    if sub.empty:
        return sub

    threshold = sub["Max_SS"].quantile(ss_quantile)
    sub = sub[sub["Max_SS"] >= threshold]
    if sub.empty:
        return sub

    sub = sub.copy()
    sub["_hub"] = (
        sub["Best_Sink_X"].round(4).astype(str) + "_" + sub["Best_Sink_Y"].round(4).astype(str)
    )
    hub_order = sub["_hub"].value_counts().head(max_hubs).index.tolist()
    parts: list[pd.DataFrame] = []
    for hi, hub in enumerate(hub_order):
        chunk = sub[sub["_hub"] == hub].nlargest(max_per_hub, "Max_SS").copy()
        chunk["hub_idx"] = hi
        chunk["line_color"] = HUB_LINE_COLORS[hi % len(HUB_LINE_COLORS)]
        parts.append(chunk)
    if not parts:
        return sub.iloc[0:0]
    out = pd.concat(parts, ignore_index=True)
    return out.drop(columns=["_hub"], errors="ignore")


def _trim_distant_hubs(pairs: pd.DataFrame, max_span_deg: float = 5.5) -> pd.DataFrame:
    """Keep one hub cluster if top hubs span too much of China (lines would vanish)."""
    if pairs.empty or pairs["hub_idx"].nunique() <= 1:
        return pairs
    xs = pairs["Best_Sink_X"].to_numpy()
    ys = pairs["Best_Sink_Y"].to_numpy()
    if (xs.max() - xs.min() > max_span_deg) or (ys.max() - ys.min() > max_span_deg):
        best = pairs.groupby("hub_idx").size().idxmax()
        return pairs[pairs["hub_idx"] == best].copy()
    return pairs


def _pairs_extent(pairs: pd.DataFrame, *, pad_ratio: float = 0.14, min_span: float = 0.6):
    xs = np.r_[pairs["Industry_X"], pairs["Best_Sink_X"]]
    ys = np.r_[pairs["Industry_Y"], pairs["Best_Sink_Y"]]
    xmin, xmax = float(xs.min()), float(xs.max())
    ymin, ymax = float(ys.min()), float(ys.max())
    w = max(xmax - xmin, min_span) * (1 + pad_ratio)
    h = max(ymax - ymin, min_span) * (1 + pad_ratio)
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    return (cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2)


def _draw_clipped_basemap(ax, china, provinces, extent):
    xmin, xmax, ymin, ymax = extent
    window = gpd.GeoDataFrame(geometry=[shapely_box(xmin, ymin, xmax, ymax)], crs="EPSG:4326")
    gpd.clip(china, window).boundary.plot(ax=ax, color="0.35", linewidth=0.9, zorder=1)
    gpd.clip(provinces, window).boundary.plot(ax=ax, color="0.7", linewidth=0.45, zorder=2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def _draw_hub_cluster(ax, pairs: pd.DataFrame):
    """Hub-and-spoke: straight lines from each source to its storage hub (SF-style)."""
    if pairs.empty:
        return

    for hub_idx in sorted(pairs["hub_idx"].unique()):
        sub = pairs[pairs["hub_idx"] == hub_idx]
        color = sub["line_color"].iloc[0]
        x0 = sub["Industry_X"].to_numpy()
        y0 = sub["Industry_Y"].to_numpy()
        x1 = sub["Best_Sink_X"].to_numpy()
        y1 = sub["Best_Sink_Y"].to_numpy()
        segments = np.stack(
            [np.column_stack([x0, y0]), np.column_stack([x1, y1])],
            axis=1,
        )
        ax.add_collection(
            LineCollection(
                segments, colors=color, linewidths=2.0, alpha=0.92, zorder=3, capstyle="round",
            )
        )
        ax.scatter(
            x0, y0, s=38, c=color, edgecolors="black", linewidths=0.45, zorder=4,
        )

    hubs = pairs.drop_duplicates(subset=["Best_Sink_X", "Best_Sink_Y", "hub_idx"])
    for row in hubs.itertuples():
        ax.scatter(
            row.Best_Sink_X, row.Best_Sink_Y,
            s=160, c=COLOR_HUB_MARKER, marker="^", edgecolors="black", linewidths=0.7, zorder=6,
        )


def _cluster_legend(ss_quantile: float, n_hubs: int) -> list[Line2D]:
    pct = int(round((1 - ss_quantile) * 100))
    items = [
        Line2D([0], [0], color=HUB_LINE_COLORS[0], linewidth=2.5,
               label=f"Source→sink spoke (top {pct}% SS_ij)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor=COLOR_HUB_MARKER,
               markeredgecolor="k", markersize=10, linestyle="", label="Storage hub (sink)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=HUB_LINE_COLORS[1],
               markeredgecolor="k", markersize=7, linestyle="", label="Industry source"),
    ]
    for i in range(min(n_hubs, len(HUB_LINE_COLORS))):
        items.append(
            Line2D([0], [0], color=HUB_LINE_COLORS[i], linewidth=2,
                   label=f"Hub cluster {i + 1}"),
        )
    return items


def plot_cluster_hub_panels(
    years: list[int],
    categories: list[str],
    output_path: str,
    *,
    ss_quantile: float = 0.85,
    min_link_km: float = 8.0,
    max_hubs: int = 3,
    max_per_hub: int = 35,
    pool_size: int = 8_000,
    admin0_path=NE_ADMIN0,
    admin1_path=NE_ADMIN1,
    dpi: int = 300,
) -> None:
    china, provinces, china_extent = _china_layers(admin0_path, admin1_path)

    ncols = 3
    nrows = math.ceil(len(years) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.8 * ncols, 6.5 * nrows))
    mf = manuscript_font_bundle(fig.get_figwidth())
    axes = np.atleast_1d(axes).flatten()
    letters = "abcdefghijklmnopqrstuvwxyz"
    max_hub_seen = 0

    for i, year in enumerate(years):
        ax = axes[i]
        try:
            pool, region_name = _load_best_region_pool(
                year, categories, min_link_km=min_link_km, pool_size=pool_size,
            )
            pairs = _select_hub_spokes(
                pool,
                ss_quantile=ss_quantile,
                min_link_km=min_link_km,
                max_hubs=max_hubs,
                max_per_hub=max_per_hub,
            )
            pairs = _trim_distant_hubs(pairs)
            if pairs.empty:
                _draw_clipped_basemap(ax, china, provinces, china_extent)
                note = f"{region_name} · no spokes"
            else:
                extent = _pairs_extent(pairs)
                _draw_clipped_basemap(ax, china, provinces, extent)
                _draw_hub_cluster(ax, pairs)
                n_hubs = pairs["hub_idx"].nunique()
                max_hub_seen = max(max_hub_seen, n_hubs)
                med = pairs["Best_Distance_km"].median()
                note = f"{region_name} · {len(pairs)} spokes · med {med:.0f} km"
        except FileNotFoundError:
            _draw_clipped_basemap(ax, china, provinces, china_extent)
            note = "no score data"

        ax.text(
            0.03, 0.97, f"({letters[i]}) {year}\n{note}",
            transform=ax.transAxes, va="top", ha="left", fontsize=mf["title"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.92, edgecolor="0.7"),
            zorder=15,
        )

    for j in range(len(years), len(axes)):
        axes[j].axis("off")

    fig.legend(
        handles=_cluster_legend(ss_quantile, max_hub_seen),
        loc="lower center", ncol=2, fontsize=mf["legend"], frameon=True, bbox_to_anchor=(0.5, 0.01),
    )
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] Cluster hub map saved to {output_path}")


def _flow_segments(matches: pd.DataFrame) -> np.ndarray:
    x0 = matches["Industry_X"].to_numpy(dtype=np.float64)
    y0 = matches["Industry_Y"].to_numpy(dtype=np.float64)
    x1 = matches["Best_Sink_X"].to_numpy(dtype=np.float64)
    y1 = matches["Best_Sink_Y"].to_numpy(dtype=np.float64)
    return np.stack(
        [np.column_stack([x0, y0]), np.column_stack([x1, y1])],
        axis=1,
    )


def _china_layers(admin0_path, admin1_path):
    admin0 = gpd.read_file(admin0_path).to_crs("EPSG:4326")
    admin1 = gpd.read_file(admin1_path).to_crs("EPSG:4326")

    china0 = admin0[admin0["ADMIN"] == "China"]
    taiwan0 = admin0[admin0["ADMIN"].astype(str).str.contains("Taiwan", case=False, na=False)]
    china = gpd.GeoDataFrame(
        pd.concat([china0, taiwan0] if not taiwan0.empty else [china0], ignore_index=True),
        crs="EPSG:4326",
    )
    if "admin" in admin1.columns:
        provinces = admin1[admin1["admin"] == "China"]
    else:
        provinces = admin1[admin1["ADM0NAME"] == "China"]

    minx, miny, maxx, maxy = china.total_bounds
    pad_x = (maxx - minx) * 0.02
    pad_y = (maxy - miny) * 0.02
    extent = (minx - pad_x, maxx + pad_x, miny - pad_y, maxy + pad_y)
    return china, provinces, extent


def _draw_china_basemap(ax, china, provinces, extent):
    xmin, xmax, ymin, ymax = extent
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    china.boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
    provinces.boundary.plot(ax=ax, color="0.65", linewidth=0.4, zorder=2)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def _draw_national_panel(
    ax,
    matches: pd.DataFrame,
    industry_bg: pd.DataFrame,
    *,
    min_link_km: float = 25.0,
):
    """Full-China panel: context industries + source→sink flow lines."""
    if not industry_bg.empty:
        ax.scatter(
            industry_bg["Industry_X"], industry_bg["Industry_Y"],
            s=3, c=COLOR_INDUSTRY_BG, alpha=0.18, linewidths=0, zorder=2, rasterized=True,
        )

    if matches.empty:
        return

    dist_km = matches["Best_Distance_km"].to_numpy(dtype=np.float64)
    has_line = dist_km >= min_link_km
    line_df = matches[has_line]
    coloc_df = matches[~has_line]

    sinks = matches.drop_duplicates(subset=["Best_Sink_X", "Best_Sink_Y"])
    dsa = sinks[sinks["Sink_StorageType"] == "DSA"]
    eor = sinks[sinks["Sink_StorageType"] == "EOR"]
    if not dsa.empty:
        ax.scatter(
            dsa["Best_Sink_X"], dsa["Best_Sink_Y"],
            s=12, c=COLOR_DSA, marker="h", edgecolors="black", linewidths=0.2, zorder=4,
        )
    if not eor.empty:
        ax.scatter(
            eor["Best_Sink_X"], eor["Best_Sink_Y"],
            s=14, c=COLOR_EOR, marker="H", edgecolors="black", linewidths=0.2, zorder=4,
        )

    if not coloc_df.empty:
        ax.scatter(
            coloc_df["Industry_X"], coloc_df["Industry_Y"],
            s=10, c=COLOR_INDUSTRY_MATCH, edgecolors="none", alpha=0.45, zorder=5, rasterized=True,
        )

    if not line_df.empty:
        ax.scatter(
            line_df["Industry_X"], line_df["Industry_Y"],
            s=14, c=COLOR_INDUSTRY_MATCH, edgecolors="white", linewidths=0.25,
            alpha=0.95, zorder=8, rasterized=True,
        )
        segments = _flow_segments(line_df)
        ss = line_df["Max_SS"].to_numpy()
        ss_norm = (ss - ss.min()) / (ss.max() - ss.min() + 1e-12)
        dist_norm = line_df["Best_Distance_km"].to_numpy()
        dist_norm = (dist_norm - dist_norm.min()) / (dist_norm.max() - dist_norm.min() + 1e-12)
        lw = 2.0 + 3.0 * ss_norm + 2.5 * dist_norm

        ax.add_collection(
            LineCollection(
                segments, colors=COLOR_FLOW_HALO, linewidths=lw + 1.8,
                alpha=0.9, zorder=9, capstyle="round",
            )
        )
        ax.add_collection(
            LineCollection(
                segments, colors=COLOR_FLOW, linewidths=lw,
                alpha=1.0, zorder=10, capstyle="round",
            )
        )
        ax.scatter(
            line_df["Best_Sink_X"], line_df["Best_Sink_Y"],
            s=16,
            c=np.where(line_df["Sink_StorageType"].to_numpy() == "EOR", COLOR_EOR, COLOR_DSA),
            marker="h", edgecolors="black", linewidths=0.25, zorder=11,
        )


def _national_legend(
    categories: list[str], ss_quantile: float, min_link_km: float,
) -> list[Line2D]:
    pct = int(round((1 - ss_quantile) * 100))
    return [
        Line2D([0], [0], color=COLOR_FLOW, linewidth=2.5,
               label=f"Source→sink flow (top {pct}% SS_ij, ≥{min_link_km:g} km)"),
        Line2D([0], [0], marker="h", color="w", markerfacecolor=COLOR_DSA,
               markeredgecolor="k", markersize=7, label="DSA storage"),
        Line2D([0], [0], marker="H", color="w", markerfacecolor=COLOR_EOR,
               markeredgecolor="k", markersize=7, label="EOR storage"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_INDUSTRY_MATCH,
               markeredgecolor="none", markersize=5, linestyle="",
               label="Industry (co-located match)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_INDUSTRY_BG,
               markeredgecolor="none", markersize=5, alpha=0.7, linestyle="",
               label="Industry (all scored, sample)"),
    ] + [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_INDUSTRY_MATCH,
               markeredgecolor="none", markersize=5, linestyle="", label=cat)
        for cat in sorted(set(categories))[:4]
    ]


def plot_national_map_panels(
    years: list[int],
    categories: list[str],
    output_path: str,
    *,
    ss_quantile: float = 0.85,
    max_flows: int = 2500,
    min_link_km: float = 25.0,
    industry_sample: int = 12_000,
    admin0_path=NE_ADMIN0,
    admin1_path=NE_ADMIN1,
    dpi: int = 300,
) -> None:
    china, provinces, extent = _china_layers(admin0_path, admin1_path)

    ncols = 3
    nrows = math.ceil(len(years) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 6.8 * nrows))
    mf = manuscript_font_bundle(fig.get_figwidth())
    axes = np.atleast_1d(axes).flatten()
    letters = "abcdefghijklmnopqrstuvwxyz"
    cats_seen: list[str] = []

    for i, year in enumerate(years):
        ax = axes[i]
        _draw_china_basemap(ax, china, provinces, extent)
        try:
            raw = _thin_matches_from_files(
                year, categories, ss_quantile, max_flows, min_link_km=min_link_km,
            )
            matches = _enrich_sink_coords(raw, year)
            bg = _sample_industries(year, categories, industry_sample, seed=year)
            _draw_national_panel(ax, matches, bg, min_link_km=min_link_km)
            cats_seen.extend(matches["Category"].tolist())
            if not matches.empty:
                flows = matches[matches["Best_Distance_km"] >= min_link_km]
                med = flows["Best_Distance_km"].median() if not flows.empty else 0.0
                note = (
                    f"{len(flows):,} flows · median {med:.0f} km · "
                    f"{len(bg):,} industries (sample)"
                )
            else:
                note = f"0 flows · {len(bg):,} industries (sample)"
        except FileNotFoundError:
            note = "no score data"

        ax.text(
            0.02, 0.98, f"({letters[i]}) {year}\n{note}",
            transform=ax.transAxes, va="top", ha="left", fontsize=mf["title"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.9, edgecolor="0.7"),
            zorder=10,
        )

    for j in range(len(years), len(axes)):
        axes[j].axis("off")

    fig.legend(
        handles=_national_legend(cats_seen, ss_quantile, min_link_km),
        loc="lower center", ncol=3, fontsize=mf["legend"], frameon=True, bbox_to_anchor=(0.5, 0.01),
    )
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] Cluster map saved to {output_path}")


# ---- Regional flow maps (optional) ----

def _load_regional_pairs(year: int, categories: list[str], box: tuple[float, float, float, float]) -> pd.DataFrame:
    sink_cols = ["Best_Sink_X", "Best_Sink_Y", "Sink_StorageType"]
    parts: list[pd.DataFrame] = []
    for short in categories:
        p = SS_SCORES_DIR / f"ss_{short}_{year}.csv"
        if not p.exists():
            continue
        header = pd.read_csv(p, nrows=0).columns.tolist()
        has_sink = all(c in header for c in sink_cols)
        usecols = [
            "Industry_X", "Industry_Y", "Max_SS", "Category", "Best_Distance_km",
            *(sink_cols if has_sink else []),
        ]
        usecols = [c for c in usecols if c in header]
        for chunk in pd.read_csv(p, usecols=usecols, chunksize=300_000):
            chunk = chunk[chunk["Max_SS"] > 0]
            chunk = chunk[_in_box(chunk, box)]
            if not chunk.empty:
                parts.append(chunk)
    if not parts:
        raise FileNotFoundError(f"No scores in region for year={year}.")
    merged = pd.concat(parts, ignore_index=True)
    merged = _enrich_sink_coords(merged, year)
    merged["SS_ij"] = merged["Max_SS"]
    merged["_cell"] = merged["Industry_X"].round(4).astype(str) + "_" + merged["Industry_Y"].round(4).astype(str)
    return merged.sort_values("SS_ij").drop_duplicates("_cell", keep="last").drop(columns="_cell")


def _filter_high_priority(
    df: pd.DataFrame, *, ss_quantile: float, min_link_km: float, max_links: int | None,
) -> pd.DataFrame:
    sub = df[(df["Industry_X"] != df["Best_Sink_X"]) | (df["Industry_Y"] != df["Best_Sink_Y"])].copy()
    if "Best_Distance_km" in sub.columns:
        sub = sub[sub["Best_Distance_km"] >= min_link_km]
    if sub.empty:
        return sub
    threshold = sub["SS_ij"].quantile(ss_quantile)
    sub = sub[sub["SS_ij"] >= threshold].sort_values("SS_ij", ascending=False)
    if max_links and len(sub) > max_links:
        sub = sub.head(max_links)
    return sub


def plot_regional_flow_map(
    year: int,
    categories: list[str],
    region: str,
    box: tuple[float, float, float, float],
    output_path: str,
    *,
    ss_quantile: float = 0.85,
    min_link_km: float = 10.0,
    max_links: int = 80,
    admin1_path=NE_ADMIN1,
    dpi: int = 300,
) -> None:
    provinces = gpd.read_file(admin1_path).to_crs("EPSG:4326")
    if "admin" in provinces.columns:
        provinces = provinces[provinces["admin"] == "China"]
    else:
        provinces = provinces[provinces["ADM0NAME"] == "China"]

    raw = _load_regional_pairs(year, categories, box)
    pairs = _filter_high_priority(raw, ss_quantile=ss_quantile, min_link_km=min_link_km, max_links=max_links)
    lines = gpd.GeoDataFrame(
        pairs,
        geometry=[
            LineString([(r.Industry_X, r.Industry_Y), (r.Best_Sink_X, r.Best_Sink_Y)])
            for r in pairs.itertuples()
        ],
        crs="EPSG:4326",
    ).to_crs(CHINA_ALBERS)
    sources = gpd.GeoDataFrame(
        pairs.drop_duplicates(subset=["Industry_X", "Industry_Y"]),
        geometry=[Point(r.Industry_X, r.Industry_Y) for r in pairs.drop_duplicates(subset=["Industry_X", "Industry_Y"]).itertuples()],
        crs="EPSG:4326",
    ).to_crs(CHINA_ALBERS)
    sinks = gpd.GeoDataFrame(
        pairs.drop_duplicates(subset=["Best_Sink_X", "Best_Sink_Y"]),
        geometry=[Point(r.Best_Sink_X, r.Best_Sink_Y) for r in pairs.drop_duplicates(subset=["Best_Sink_X", "Best_Sink_Y"]).itertuples()],
        crs="EPSG:4326",
    ).to_crs(CHINA_ALBERS)
    prov = provinces.to_crs(CHINA_ALBERS)

    fig, ax = plt.subplots(figsize=(8, 7))
    mf = manuscript_font_bundle(fig.get_figwidth())
    bounds = lines.total_bounds if not lines.empty else box
    if not lines.empty:
        xmin, ymin, xmax, ymax = lines.total_bounds
        pad_x, pad_y = (xmax - xmin) * 0.1, (ymax - ymin) * 0.1
        extent = (xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y)
    else:
        extent = prov.total_bounds

    window = gpd.GeoDataFrame(geometry=[shapely_box(*extent)], crs=CHINA_ALBERS)
    gpd.clip(prov, window).boundary.plot(ax=ax, color="0.55", linewidth=0.6, zorder=1)
    if not lines.empty:
        lines.plot(ax=ax, color=COLOR_FLOW_REGIONAL, linewidth=1.3, zorder=3)
        sources.plot(ax=ax, color=COLOR_SOURCE_REGIONAL, markersize=50, edgecolor="k", zorder=5)
        dsa = sinks[sinks["Sink_StorageType"] == "DSA"]
        eor = sinks[sinks["Sink_StorageType"] == "EOR"]
        if not dsa.empty:
            dsa.plot(ax=ax, color=COLOR_DSA, marker="p", markersize=80, edgecolor="k", zorder=6)
        if not eor.empty:
            eor.plot(ax=ax, color=COLOR_EOR, marker="*", markersize=90, edgecolor="k", zorder=6)

    ax.set_xlim(extent[0], extent[2])
    ax.set_ylim(extent[1], extent[3])
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(0.03, 0.97, f"{year} · {region}\n{len(pairs)} flows",
            transform=ax.transAxes, va="top", fontsize=mf["title"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.9))
    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] Regional flow map saved to {output_path}")


def parse_args():
    p = argparse.ArgumentParser(description="China source–sink cluster map.")
    p.add_argument(
        "--layout",
        choices=("national", "hub", "regional"),
        default="national",
        help="national=full China (default); hub=zoomed spokes; regional=one macro-region.",
    )
    p.add_argument("--years", type=int, nargs="+", default=DEFAULT_MAP_YEARS,
                   help=f"Panel years (default: {DEFAULT_MAP_YEARS}).")
    p.add_argument("--year", type=int, default=None)
    p.add_argument("--categories", nargs="+", default=["MC"])
    p.add_argument("--ss-quantile", type=float, default=0.85)
    p.add_argument("--max-flows", type=int, default=2500)
    p.add_argument("--industry-sample", type=int, default=12_000)
    p.add_argument(
        "--min-link-km", type=float, default=25.0,
        help="Minimum source→sink distance (km) for drawn flows; default 25 for national.",
    )
    p.add_argument("--max-hubs", type=int, default=3)
    p.add_argument("--max-per-hub", type=int, default=35)
    p.add_argument("--region", default="West", choices=list(MATCH_REGIONS.keys()))
    p.add_argument("--output", default=None)
    return p.parse_args()


def main():
    args = parse_args()
    ensure_charts_dir()
    cats = _resolve_categories(args.categories)
    years = [args.year] if args.year is not None else args.years

    if args.layout == "regional":
        year = args.year or args.years[0]
        out = args.output or chart_path(f"source_sink_flows_{year}_{args.region}.png")
        plot_regional_flow_map(
            year, cats, args.region, MATCH_REGIONS[args.region], out,
            ss_quantile=args.ss_quantile, min_link_km=args.min_link_km,
            max_links=min(args.max_flows, 80),
        )
        return

    if args.layout == "hub":
        out = args.output or chart_path("source_sink_cluster_hubs.png")
        plot_cluster_hub_panels(
            years, cats, out,
            ss_quantile=args.ss_quantile,
            min_link_km=args.min_link_km,
            max_hubs=args.max_hubs,
            max_per_hub=args.max_per_hub,
        )
        return

    out = args.output or chart_path("source_sink_cluster_map.png")
    plot_national_map_panels(
        years, cats, out,
        ss_quantile=args.ss_quantile,
        max_flows=args.max_flows,
        min_link_km=args.min_link_km,
        industry_sample=args.industry_sample,
    )


if __name__ == "__main__":
    main()
