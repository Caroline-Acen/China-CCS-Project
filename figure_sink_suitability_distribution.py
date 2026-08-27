"""
Sink suitability score (SS) distribution figure — alternative to storage_industry_map.

Method (Source-Sink Matching.docx / DIDPA):
    SS_{i,j} = (R_{i,t} * S_{i,t}) / d_{i,j}^2
    R, S in tonnes and tonnes/year; d in km (great-circle).
    For each industry grid j, consider storage cells i within 250 km and
    record the maximum SS (best sink for that grid).

Workflow (strictly sequential, one industry group × year at a time):
    Sinks: all four data grid workbooks — DSA/EOR × Pipeline/Tank (X, Y, Storage
    Potential, Injection Rate; R_yyyy / S_yyyy for later years). Subsurface R and S
    define SS; transport mode is included so both delivery options are in the sink pool.
    1. Compute → data/ss_scores/ss_{TC|ME|…}_{year}.csv
    2. Merge   → data/ss_scores/sink_suitability_combined.csv
    3. Plot    → charts/sink_suitability_pareto.png (or _bivariate.png)
               → data/sink_suitability_summary.csv

Does not modify figure4a.py / storage_industry_map.png.
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.neighbors import BallTree

from config import (
    COST_YEARS,
    NE_ADMIN0,
    NE_ADMIN1,
    INDUSTRIES_CSV,
    SS_SCORES_DIR,
    DSA_PIPELINE,
    DSA_TANK,
    EOR_PIPELINE,
    EOR_TANK,
    chart_path,
    data_path,
    ensure_charts_dir,
    manuscript_font_bundle,
)
from didpa import (
    EARTH_RADIUS_KM,
    MT_TO_T,
    load_storage_sinks,
    sink_suitability_score,
)

MAX_PAIR_DISTANCE_KM = 250.0
COMBINED_CSV_NAME = "sink_suitability_combined.csv"

INDUSTRY_CATEGORIES = {
    "TC19": "TC",
    "ME19": "ME",
    "WF19": "WF",
    "AF19": "AF",
    "OM19": "OM",
    "MC19": "MC",
    "PM19": "PM",
    "PP19": "PP",
}

CATEGORY_KEYS = list(INDUSTRY_CATEGORIES.keys())


def _log(msg: str, t0: float) -> None:
    elapsed = time.perf_counter() - t0
    print(f"[{elapsed:8.1f}s] {msg}")


def _ensure_ss_dir() -> Path:
    SS_SCORES_DIR.mkdir(parents=True, exist_ok=True)
    return SS_SCORES_DIR


def _category_score_path(category_col: str, year: int) -> Path:
    short = INDUSTRY_CATEGORIES[category_col]
    return _ensure_ss_dir() / f"ss_{short}_{year}.csv"


def combined_csv_path() -> Path:
    return _ensure_ss_dir() / COMBINED_CSV_NAME


def _sink_sources(storage_mode: str) -> list[tuple[Path, str, str]]:
    """
    Return (excel_path, storage_type, transport_mode) for the requested pool.

    storage_mode:
        all      — DSA/EOR pipeline + tank (default)
        pipeline — DSA + EOR pipeline only
        tank     — DSA + EOR tank only
        dsa      — DSA pipeline + tank
        eor      — EOR pipeline + tank
        both     — alias for pipeline (legacy)
    """
    catalog = {
        ("DSA", "pipeline"): DSA_PIPELINE,
        ("DSA", "tank"): DSA_TANK,
        ("EOR", "pipeline"): EOR_PIPELINE,
        ("EOR", "tank"): EOR_TANK,
    }
    mode = storage_mode.lower()
    if mode in ("all", "full"):
        keys = list(catalog.keys())
    elif mode == "pipeline":
        keys = [("DSA", "pipeline"), ("EOR", "pipeline")]
    elif mode == "tank":
        keys = [("DSA", "tank"), ("EOR", "tank")]
    elif mode == "dsa":
        keys = [("DSA", "pipeline"), ("DSA", "tank")]
    elif mode in ("eor",):
        keys = [("EOR", "pipeline"), ("EOR", "tank")]
    elif mode == "both":
        keys = [("DSA", "pipeline"), ("EOR", "pipeline")]
    else:
        raise ValueError(
            f"Unknown storage_mode {storage_mode!r}. "
            "Use all, pipeline, tank, dsa, eor, or both."
        )
    return [(catalog[k], k[0], k[1]) for k in keys]


def load_sinks_for_year(
    year: int,
    storage_mode: str = "all",
    admin1_dir=None,
) -> pd.DataFrame:
    """Load storage grids (pipeline and/or tank); R,S in tonnes for SS."""
    from config import NE_ADMIN1

    admin1_dir = admin1_dir or NE_ADMIN1
    frames = []
    for path, stype, transport in _sink_sources(storage_mode):
        sinks = load_storage_sinks(path, year, admin1_dir=admin1_dir)
        sinks["StorageType"] = stype
        sinks["TransportMode"] = transport
        sinks["SourceFile"] = path.name
        frames.append(sinks)
    if not frames:
        raise ValueError(f"No sinks loaded for storage_mode={storage_mode!r}")
    sinks = pd.concat(frames, ignore_index=True)
    # Same X,Y can appear in pipeline and tank; keep the stronger R×S cell
    sinks["_rs"] = sinks["R_mt"] * sinks["S_mt"]
    sinks = (
        sinks.sort_values("_rs", ascending=False)
        .drop_duplicates(subset=["X", "Y"], keep="first")
        .drop(columns=["_rs"])
        .reset_index(drop=True)
    )
    sinks["R_t"] = sinks["R_mt"].astype(float) * MT_TO_T
    sinks["S_t"] = sinks["S_mt"].astype(float) * MT_TO_T
    return sinks


def _max_ss_for_emitters(
    lat: np.ndarray,
    lon: np.ndarray,
    tree: BallTree,
    r_t: np.ndarray,
    s_t: np.ndarray,
    max_km: float = MAX_PAIR_DISTANCE_KM,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each emitter, max SS to any sink within max_km. Returns max_ss, best_dist_km, best_sink_idx."""
    emitters_rad = np.column_stack([np.radians(lat), np.radians(lon)])
    radius_rad = max_km / EARTH_RADIUS_KM
    neighbor_lists, dist_lists = tree.query_radius(
        emitters_rad, r=radius_rad, return_distance=True
    )

    n = len(lat)
    max_ss = np.zeros(n, dtype=np.float64)
    best_dist = np.full(n, np.nan, dtype=np.float64)
    best_sink_idx = np.full(n, -1, dtype=np.int32)

    for idx in range(n):
        sink_idx = neighbor_lists[idx]
        if len(sink_idx) == 0:
            continue
        dist_km = dist_lists[idx] * EARTH_RADIUS_KM
        best = -1.0
        best_d = np.nan
        best_si = -1
        for j, sidx in enumerate(sink_idx):
            d = float(dist_km[j])
            if d <= 0:
                continue
            ss = sink_suitability_score(float(r_t[sidx]), float(s_t[sidx]), d)
            if ss > best:
                best = ss
                best_d = d
                best_si = int(sidx)
        if best > 0:
            max_ss[idx] = best
            best_dist[idx] = best_d
            best_sink_idx[idx] = best_si

    return max_ss, best_dist, best_sink_idx


def compute_category_year(
    category_col: str,
    year: int,
    *,
    storage_mode: str = "all",
    chunk_size: int = 250_000,
    industries_csv: Path = INDUSTRIES_CSV,
    force: bool = False,
    max_chunks: int | None = None,
) -> Path:
    """SS for one industry category and one year; industries read in chunks."""
    out_path = _category_score_path(category_col, year)
    if out_path.exists() and not force:
        print(f"[skip] {out_path.name} already exists")
        return out_path

    label = INDUSTRY_CATEGORIES[category_col]
    t0 = time.perf_counter()
    print(f"\n=== {label} / {year} ===")

    sinks = load_sinks_for_year(year, storage_mode=storage_mode)
    sources = _sink_sources(storage_mode)
    _log(
        f"{len(sinks):,} storage cells for {year} "
        f"({storage_mode}: {', '.join(p.name for p, _, _ in sources)})",
        t0,
    )

    sink_x = sinks["X"].to_numpy(dtype=np.float64)
    sink_y = sinks["Y"].to_numpy(dtype=np.float64)
    sink_storage_type = sinks["StorageType"].to_numpy()
    r_t = sinks["R_t"].to_numpy()
    s_t = sinks["S_t"].to_numpy()
    tree = BallTree(
        np.column_stack([np.radians(sink_y), np.radians(sink_x)]),
        metric="haversine",
    )

    del sinks

    header_written = False
    total_rows = 0

    for chunk_idx, chunk in enumerate(
        pd.read_csv(
            industries_csv,
            usecols=["X", "Y", category_col],
            chunksize=chunk_size,
            encoding="utf-8",
            low_memory=False,
        ),
        start=1,
    ):
        chunk = chunk[chunk[category_col] > 0]
        if chunk.empty:
            if max_chunks is not None and chunk_idx >= max_chunks:
                break
            continue

        lat = chunk["Y"].to_numpy(dtype=np.float64)
        lon = chunk["X"].to_numpy(dtype=np.float64)
        max_ss, best_dist, best_si = _max_ss_for_emitters(lat, lon, tree, r_t, s_t)
        valid = max_ss > 0
        if valid.any():
            si = best_si[valid]
            out = pd.DataFrame(
                {
                    "Year": year,
                    "Category": label,
                    "CategoryColumn": category_col,
                    "Industry_X": chunk["X"].to_numpy()[valid],
                    "Industry_Y": chunk["Y"].to_numpy()[valid],
                    "Max_SS": max_ss[valid],
                    "Best_Distance_km": best_dist[valid],
                    "Best_Sink_X": sink_x[si],
                    "Best_Sink_Y": sink_y[si],
                    "Sink_StorageType": sink_storage_type[si],
                }
            )
            out.to_csv(
                out_path,
                mode="a" if header_written else "w",
                index=False,
                header=not header_written,
            )
            header_written = True
            total_rows += len(out)

        if chunk_idx % 20 == 0:
            _log(f"  industries chunk {chunk_idx}: {total_rows:,} scores saved", t0)
        if max_chunks is not None and chunk_idx >= max_chunks:
            break

    if not header_written:
        pd.DataFrame(
            columns=[
                "Year", "Category", "CategoryColumn",
                "Industry_X", "Industry_Y", "Max_SS", "Best_Distance_km",
                "Best_Sink_X", "Best_Sink_Y", "Sink_StorageType",
            ]
        ).to_csv(out_path, index=False)

    _log(f"[OK] {out_path.name} ({total_rows:,} rows)", t0)
    return out_path


def compute_all_sequential(
    years: list[int],
    categories: list[str] | None = None,
    storage_mode: str = "all",
    chunk_size: int = 250_000,
    force: bool = False,
    max_chunks: int | None = None,
) -> None:
    """One industry group and one year at a time (nested loops, no parallelism)."""
    categories = categories or CATEGORY_KEYS
    jobs = [(y, c) for y in years for c in categories]
    total_jobs = len(jobs)

    print(f"Sequential SS compute: {total_jobs} jobs "
          f"({len(categories)} categories × {len(years)} years)")

    for job_num, (year, cat) in enumerate(jobs, start=1):
        print(f"\n--- Job {job_num}/{total_jobs} ---")
        compute_category_year(
            cat,
            year,
            storage_mode=storage_mode,
            chunk_size=chunk_size,
            force=force,
            max_chunks=max_chunks,
        )


def _read_category_csv(path: Path) -> pd.DataFrame:
    if path.suffix == ".gz" or path.name.endswith(".csv.gz"):
        return pd.read_csv(path, compression="gzip")
    return pd.read_csv(path)


def merge_category_csvs(
    years: list[int],
    categories: list[str] | None = None,
    output_path: Path | None = None,
) -> Path:
    """Concatenate all per-category/year CSVs into one combined CSV."""
    categories = categories or CATEGORY_KEYS
    output_path = output_path or combined_csv_path()

    parts: list[pd.DataFrame] = []
    missing: list[str] = []

    for year in years:
        for cat in categories:
            path = _category_score_path(cat, year)
            gz_legacy = path.parent / f"{path.stem}.csv.gz"
            if path.exists():
                parts.append(_read_category_csv(path))
            elif gz_legacy.exists():
                parts.append(_read_category_csv(gz_legacy))
            else:
                missing.append(path.name)

    if missing:
        raise FileNotFoundError(
            "Missing score CSVs (run compute first):\n  " + "\n  ".join(missing[:20])
            + (f"\n  ... and {len(missing) - 20} more" if len(missing) > 20 else "")
        )

    combined = pd.concat(parts, ignore_index=True)
    combined.to_csv(output_path, index=False)
    print(f"[OK] Merged {len(parts)} files -> {output_path} ({len(combined):,} rows)")
    return output_path


def load_combined_scores(path: Path | None = None) -> pd.DataFrame:
    path = path or combined_csv_path()
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run with --merge after --compute."
        )
    return pd.read_csv(path)


def _load_scores_for_year(
    year: int,
    categories: list[str] | None = None,
) -> pd.DataFrame:
    """Load Max_SS and Best_Distance_km for one year from per-category CSVs."""
    categories = categories or CATEGORY_KEYS
    usecols = ["Max_SS", "Best_Distance_km"]
    parts: list[pd.DataFrame] = []
    for cat in categories:
        path = _category_score_path(cat, year)
        if not path.exists():
            continue
        parts.append(pd.read_csv(path, usecols=usecols))
    if not parts:
        return pd.DataFrame(columns=usecols)
    return pd.concat(parts, ignore_index=True)


def _scores_for_plotting(
    year: int,
    categories: list[str] | None,
    combined_path: Path | None,
) -> pd.DataFrame:
    if combined_path is not None and combined_path.exists():
        df = load_combined_scores(combined_path)
        return df[df["Year"] == year]
    return _load_scores_for_year(year, categories)


def _load_scores_for_year_map(
    year: int,
    categories: list[str] | None = None,
) -> pd.DataFrame:
    """Load map-ready columns for one year from per-category CSVs."""
    categories = categories or CATEGORY_KEYS
    usecols = ["Industry_X", "Industry_Y", "Max_SS", "Best_Distance_km"]
    parts: list[pd.DataFrame] = []
    for cat in categories:
        path = _category_score_path(cat, year)
        if not path.exists():
            continue
        parts.append(pd.read_csv(path, usecols=usecols))
    if not parts:
        return pd.DataFrame(columns=usecols)
    return pd.concat(parts, ignore_index=True)


def _load_year_from_combined(
    combined_path: Path,
    year: int,
    usecols: list[str],
) -> pd.DataFrame:
    """Read a single year from a large combined CSV without loading all rows."""
    parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(combined_path, usecols=["Year", *usecols], chunksize=1_000_000):
        y = chunk[chunk["Year"] == year]
        if not y.empty:
            parts.append(y[usecols])
    if not parts:
        return pd.DataFrame(columns=usecols)
    return pd.concat(parts, ignore_index=True)


def _scores_for_map(
    year: int,
    categories: list[str] | None,
    combined_path: Path | None,
) -> pd.DataFrame:
    usecols = ["Industry_X", "Industry_Y", "Max_SS", "Best_Distance_km"]
    if combined_path is not None and combined_path.exists():
        return _load_year_from_combined(combined_path, year, usecols)
    return _load_scores_for_year_map(year, categories)


def _china_layers(admin0_path=NE_ADMIN0, admin1_path=NE_ADMIN1):
    admin0 = gpd.read_file(admin0_path)
    admin1 = gpd.read_file(admin1_path)

    china = admin0[admin0["ADMIN"] == "China"].copy()
    if "admin" in admin1.columns:
        provinces = admin1[admin1["admin"] == "China"].copy()
    elif "ADM0NAME" in admin1.columns:
        provinces = admin1[admin1["ADM0NAME"] == "China"].copy()
    else:
        provinces = gpd.GeoDataFrame(geometry=[], crs=admin1.crs)

    if china.crs is not None and provinces.crs is not None and china.crs != provinces.crs:
        provinces = provinces.to_crs(china.crs)

    # Match cost-map geometry by plotting in Web Mercator.
    china = china.to_crs(epsg=3857)
    if len(provinces) > 0:
        provinces = provinces.to_crs(epsg=3857)

    xmin, ymin, xmax, ymax = china.total_bounds
    pad_x = 80_000
    pad_y = 80_000
    extent = (xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y)
    return china, provinces, extent


def _draw_map_basemap(ax, china, provinces, extent):
    ax.set_facecolor("white")
    china.boundary.plot(ax=ax, color="#666666", linewidth=0.8, zorder=1)
    if len(provinces) > 0:
        provinces.boundary.plot(ax=ax, color="#d0d0d0", linewidth=0.5, zorder=1)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal", adjustable="box")


def _box_panel(ax):
    ax.set_axis_on()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_edgecolor("#4a4a4a")


def _map_color_limits(
    years: list[int],
    categories: list[str] | None,
    combined_path: Path | None,
    sample_per_year: int = 160_000,
) -> tuple[float, float]:
    logs: list[np.ndarray] = []
    for year in years:
        df = _scores_for_map(year, categories, combined_path)
        ss = df["Max_SS"].to_numpy(dtype=np.float64)
        valid = (ss > 0) & np.isfinite(ss)
        if not valid.any():
            continue
        lv = np.log10(ss[valid])
        if len(lv) > sample_per_year:
            rng = np.random.default_rng(year)
            idx = rng.choice(len(lv), size=sample_per_year, replace=False)
            lv = lv[idx]
        logs.append(lv)

    if not logs:
        return 0.0, 1.0

    all_logs = np.concatenate(logs)
    vmin = float(np.quantile(all_logs, 0.10))
    vmax = float(np.quantile(all_logs, 0.95))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin = float(np.nanmin(all_logs))
        vmax = float(np.nanmax(all_logs))
    if vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def _map_color_limits_for_years(
    years: list[int],
    categories: list[str] | None,
    combined_path: Path | None,
) -> tuple[float, float]:
    """Compute robust log10(Max_SS) limits for a specific subset of years."""
    logs: list[np.ndarray] = []
    for year in years:
        df = _scores_for_map(year, categories, combined_path)
        ss = df["Max_SS"].to_numpy(dtype=np.float64)
        valid = (ss > 0) & np.isfinite(ss)
        if valid.any():
            logs.append(np.log10(ss[valid]))

    if not logs:
        return 0.0, 1.0

    all_logs = np.concatenate(logs)
    vmin = float(np.quantile(all_logs, 0.10))
    vmax = float(np.quantile(all_logs, 0.95))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin = float(np.nanmin(all_logs))
        vmax = float(np.nanmax(all_logs))
    if vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def plot_map_panels(
    years: list[int],
    output_path: str,
    *,
    combined_path: Path | None = None,
    categories: list[str] | None = None,
    gridsize: int = 72,
    mincnt: int = 8,
    cmap: str = "YlGn",
) -> None:
    """Hex-cell map panels: source-grid sink suitability over China (one panel per year)."""
    n_years = len(years)
    ncols = 3
    nrows = math.ceil(n_years / ncols)

    china, provinces, extent = _china_layers()
    vmin, vmax = _map_color_limits(years, categories, combined_path)

    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 12))
    mf = manuscript_font_bundle(fig.get_figwidth(), target_pt=9.0)
    axes = np.atleast_1d(axes).flatten()
    mappable_row1 = None
    mappable_row2 = None
    panel_distance_ranges: dict[int, tuple[float, float]] = {}

    row1_years = years[:3]
    row2_years = years[3:]
    row1_vmin, row1_vmax = _map_color_limits_for_years(row1_years, categories, combined_path)
    row2_vmin, row2_vmax = _map_color_limits_for_years(row2_years, categories, combined_path)

    for i, year in enumerate(years):
        ax = axes[i]
        _draw_map_basemap(ax, china, provinces, extent)

        df = _scores_for_map(year, categories, combined_path)
        points = gpd.GeoDataFrame(
            df[["Industry_X", "Industry_Y"]].copy(),
            geometry=gpd.points_from_xy(df["Industry_X"], df["Industry_Y"]),
            crs="EPSG:4326",
        ).to_crs(china.crs)

        x = points.geometry.x.to_numpy(dtype=np.float64)
        y = points.geometry.y.to_numpy(dtype=np.float64)
        ss = df["Max_SS"].to_numpy(dtype=np.float64)
        mask = (
            np.isfinite(x) & np.isfinite(y) & np.isfinite(ss)
            & (ss > 0)
            & (x >= extent[0]) & (x <= extent[1])
            & (y >= extent[2]) & (y <= extent[3])
        )

        n_valid = int(mask.sum())
        if n_valid > 0:
            if i < 3:
                vmin_use, vmax_use = row1_vmin, row1_vmax
            else:
                vmin_use, vmax_use = row2_vmin, row2_vmax

            mappable = ax.hexbin(
                x[mask],
                y[mask],
                C=np.log10(ss[mask]),
                reduce_C_function=np.nanmedian,
                gridsize=gridsize,
                mincnt=mincnt,
                cmap=cmap,
                vmin=vmin_use,
                vmax=vmax_use,
                extent=extent,
                linewidths=0.16,
                edgecolors="#f3f3f3",
                zorder=2,
            )
            if i < 3:
                mappable_row1 = mappable
            else:
                mappable_row2 = mappable

        dist = df["Best_Distance_km"].to_numpy(dtype=np.float64)
        dist_mask = mask & (dist > 0) & np.isfinite(dist)
        if dist_mask.any():
            dvals = dist[dist_mask]
            dmax_raw = float(np.quantile(dvals, 0.95))
            if not np.isfinite(dmax_raw) or dmax_raw <= 0:
                dmax_raw = float(np.nanmax(dvals)) if len(dvals) else 1.0

            # Round to publication-friendly "nice" distances.
            if dmax_raw <= 120:
                step = 20.0
            elif dmax_raw <= 260:
                step = 25.0
            elif dmax_raw <= 600:
                step = 50.0
            else:
                step = 100.0

            dmin = 0.0
            dmax = float(np.ceil(dmax_raw / step) * step)
            if dmax <= dmin:
                dmax = dmin + step
            panel_distance_ranges[i] = (dmin, dmax, step)

        ax.set_title(f"{year}", fontsize=mf["title"], fontweight="bold", pad=8)
        _box_panel(ax)

    for j in range(len(years), len(axes)):
        axes[j].axis("off")

    # Match subplot geometry used by the cost maps for visual uniformity.
    fig.subplots_adjust(left=0.04, right=0.88, hspace=0.32, wspace=0.1)

    # Add one full-width distance colorbar below each subplot (outside panel boxes).
    for i in range(len(years)):
        if i not in panel_distance_ranges:
            continue
        dmin, dmax, step = panel_distance_ranges[i]
        ax = axes[i]
        pos = ax.get_position()
        bar_left = pos.x0
        bar_width = pos.width
        # Use the same thickness as the vertical colorbar width in this figure.
        bar_height = 0.015
        bar_bottom = pos.y0 - 0.026

        dax = fig.add_axes([bar_left, bar_bottom, bar_width, bar_height])
        dsm = plt.cm.ScalarMappable(cmap=f"{cmap}_r", norm=plt.Normalize(vmin=dmin, vmax=dmax))
        dsm._A = []
        dcbar = fig.colorbar(dsm, cax=dax, orientation="horizontal")
        ticks = np.arange(0.0, 251.0, 50.0)
        dcbar.set_ticks(ticks)
        dcbar.set_ticklabels([f"{t:.0f}" for t in ticks])
        dcbar.ax.tick_params(labelsize=mf["colorbar_tick"], pad=1)
        dcbar.set_label("Distance (km)", fontsize=mf["colorbar_label"], labelpad=1)

    if mappable_row1 is not None:
        cax1 = fig.add_axes([0.91, 0.53, 0.015, 0.35])
        cbar1 = fig.colorbar(mappable_row1, cax=cax1)
        cbar1.set_label("Sink suitability score", fontsize=mf["colorbar_label"])
        cbar1.ax.tick_params(labelsize=mf["colorbar_tick"])

    if mappable_row2 is not None:
        cax2 = fig.add_axes([0.91, 0.089, 0.015, 0.35])
        cbar2 = fig.colorbar(mappable_row2, cax=cax2)
        cbar2.set_label("Sink suitability score", fontsize=mf["colorbar_label"])
        cbar2.ax.tick_params(labelsize=mf["colorbar_tick"])

    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] Map plot saved to {output_path}")


def summarize_scores_by_year(
    years: list[int],
    categories: list[str] | None = None,
    combined_path: Path | None = None,
    *,
    ss_quantile: float = 0.85,
) -> pd.DataFrame:
    """Per-year summary stats for the SS figures (Pareto + bivariate)."""
    letters = "abcdefghijklmnopqrstuvwxyz"
    rows: list[dict] = []

    for i, year in enumerate(years):
        df = _scores_for_plotting(year, categories, combined_path)
        ss = df["Max_SS"].to_numpy(dtype=np.float64)
        dist = df["Best_Distance_km"].to_numpy(dtype=np.float64)
        valid = (ss > 0) & np.isfinite(ss) & np.isfinite(dist) & (dist > 0)
        ss_v = ss[valid]
        dist_v = dist[valid]

        row: dict = {
            "panel": f"({letters[i]})",
            "year": year,
            "n_grids": int(len(ss_v)),
        }
        if ss_v.size == 0:
            rows.append(row)
            continue

        q85 = float(np.quantile(ss_v, ss_quantile))
        top = ss_v >= q85
        row.update(
            {
                "max_ss_min": float(ss_v.min()),
                "max_ss_p25": float(np.quantile(ss_v, 0.25)),
                "max_ss_median": float(np.median(ss_v)),
                "max_ss_p75": float(np.quantile(ss_v, 0.75)),
                "max_ss_p85_threshold": q85,
                "max_ss_p95": float(np.quantile(ss_v, 0.95)),
                "max_ss_max": float(ss_v.max()),
                "best_distance_km_min": float(dist_v.min()),
                "best_distance_km_median": float(np.median(dist_v)),
                "best_distance_km_mean": float(dist_v.mean()),
                "best_distance_km_p85": float(np.quantile(dist_v, 0.85)),
                "best_distance_km_max": float(dist_v.max()),
                "n_top15pct_ss": int(top.sum()),
                "pct_top15pct_ss": round(100.0 * top.mean(), 2),
                "n_distance_ge_25km": int((dist_v >= 25.0).sum()),
                "pct_distance_ge_25km": round(100.0 * (dist_v >= 25.0).mean(), 2),
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def write_plot_summary_csv(
    years: list[int],
    output_path: str | Path,
    *,
    categories: list[str] | None = None,
    combined_path: Path | None = None,
) -> Path:
    summary = summarize_scores_by_year(years, categories, combined_path)
    output_path = Path(output_path)
    summary.to_csv(output_path, index=False)
    print(f"[OK] Summary CSV saved to {output_path}")
    return output_path


def _log10_axis_ticks(ax, axis: str, values: np.ndarray, *, step: float = 1.0) -> None:
    """Label a linear log10-transformed axis as powers of ten."""
    lo = float(np.floor(np.log10(values.min())))
    hi = float(np.ceil(np.log10(values.max())))
    ticks = np.arange(lo, hi + step * 0.5, step)
    labels = [
        rf"$10^{{{int(t)}}}$" if t == int(t) else rf"$10^{{{t:g}}}$"
        for t in ticks
    ]
    setter = ax.set_xticks if axis == "x" else ax.set_yticks
    labeler = ax.set_xticklabels if axis == "x" else ax.set_yticklabels
    setter(ticks)
    labeler(labels)


def plot_pareto_panels(
    years: list[int],
    output_path: str,
    *,
    combined_path: Path | None = None,
    categories: list[str] | None = None,
    log_y: bool = True,
) -> None:
    """Ordered Pareto: rank vs Max_SS (one subplot per year)."""
    n_years = len(years)
    ncols = 3
    nrows = math.ceil(n_years / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows))
    mf = manuscript_font_bundle(fig.get_figwidth())
    axes = np.atleast_1d(axes).flatten()
    letters = "abcdefghijklmnopqrstuvwxyz"

    for i, year in enumerate(years):
        ax = axes[i]
        df = _scores_for_plotting(year, categories, combined_path)
        if df.empty:
            ax.set_title(f"({letters[i]}) {year} (no data)")
            continue
        scores = np.sort(df["Max_SS"].to_numpy(dtype=np.float64))
        scores = scores[scores > 0]
        if scores.size == 0:
            ax.set_title(f"({letters[i]}) {year} (no data)")
            continue
        scores = scores[::-1]
        ranks = np.arange(1, len(scores) + 1)

        ax.plot(ranks, scores, color="#1f77b4", linewidth=1.2, alpha=0.85)
        if log_y:
            ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlabel("Rank (descending SS)", fontsize=mf["label"])
        ax.set_ylabel("Max sink suitability score", fontsize=mf["label"])
        ax.set_title(f"({letters[i]}) {year}", fontsize=mf["title"], fontweight="bold")
        ax.tick_params(axis="both", labelsize=mf["tick"])

    for j in range(len(years), len(axes)):
        axes[j].axis("off")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] Pareto plot saved to {output_path}")


def plot_bivariate_panels(
    years: list[int],
    output_path: str,
    *,
    combined_path: Path | None = None,
    categories: list[str] | None = None,
    gridsize: int = 50,
) -> None:
    """Bivariate hexbin density: best distance vs Max_SS (log scales, one panel/year)."""
    n_years = len(years)
    ncols = 3
    nrows = math.ceil(n_years / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.8 * nrows))
    mf = manuscript_font_bundle(fig.get_figwidth())
    axes = np.atleast_1d(axes).flatten()
    letters = "abcdefghijklmnopqrstuvwxyz"

    for i, year in enumerate(years):
        ax = axes[i]
        df = _scores_for_plotting(year, categories, combined_path)
        x = df["Best_Distance_km"].to_numpy(dtype=np.float64)
        y = df["Max_SS"].to_numpy(dtype=np.float64)
        mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
        if not mask.any():
            ax.set_title(f"({letters[i]}) {year} (no data)")
            continue

        x, y = x[mask], y[mask]
        lx = np.log10(x)
        ly = np.log10(y)
        hb = ax.hexbin(
            lx, ly,
            gridsize=gridsize,
            bins="log",
            cmap="viridis",
            mincnt=1,
        )
        ax.set_xlabel("Distance to best sink (km)", fontsize=mf["label"])
        ax.set_ylabel("Max sink suitability score", fontsize=mf["label"])
        ax.set_title(
            f"({letters[i]}) {year}  (n = {len(x):,})",
            fontsize=mf["title"],
            fontweight="bold",
        )
        _log10_axis_ticks(ax, "x", x, step=1.0)
        _log10_axis_ticks(ax, "y", y, step=2.0)
        ax.set_xlim(lx.min() - 0.15, min(lx.max() + 0.15, np.log10(250.0)))
        ax.set_ylim(ly.min() - 0.15, ly.max() + 0.15)
        cbar = fig.colorbar(hb, ax=ax, shrink=0.85)
        cbar.set_label("Grid count (log)", fontsize=mf["colorbar_label"])
        cbar.ax.tick_params(labelsize=mf["colorbar_tick"])
        ax.tick_params(axis="both", labelsize=mf["tick"])

    for j in range(len(years), len(axes)):
        axes[j].axis("off")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[OK] Bivariate plot saved to {output_path}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sequential SS by industry/year → CSVs → merge → plot.",
    )
    p.add_argument("--compute", action="store_true", help="Step 1: score each category/year.")
    p.add_argument("--merge", action="store_true", help="Step 2: merge CSVs into one file.")
    p.add_argument("--plot", action="store_true", help="Step 3: figure from merged CSV.")
    p.add_argument(
        "--plot-type",
        choices=("pareto", "bivariate", "map"),
        default="pareto",
    )
    p.add_argument("--years", type=int, nargs="+", default=list(COST_YEARS))
    p.add_argument("--categories", nargs="+", default=CATEGORY_KEYS)
    p.add_argument("--storage", choices=("all", "pipeline", "tank", "dsa", "eor", "both"), default="all")
    p.add_argument("--chunk-size", type=int, default=250_000)
    p.add_argument("--force", action="store_true", help="Overwrite existing category CSVs.")
    p.add_argument("--max-chunks", type=int, default=None, help="Limit chunks (testing only).")
    p.add_argument("--output-plot", default=None)
    p.add_argument(
        "--combined-csv",
        default=None,
        help="Path to merged CSV for --plot (default: data/ss_scores/sink_suitability_combined.csv).",
    )
    p.add_argument(
        "--output-summary",
        default=None,
        help="Summary CSV path (default: data/sink_suitability_summary.csv).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ensure_charts_dir()
    combined_path = Path(args.combined_csv) if args.combined_csv else None

    if not (args.compute or args.merge or args.plot):
        args.compute = True
        args.merge = True
        args.plot = True

    if args.compute:
        compute_all_sequential(
            args.years,
            categories=args.categories,
            storage_mode=args.storage,
            chunk_size=args.chunk_size,
            force=args.force,
            max_chunks=args.max_chunks,
        )

    if args.merge:
        merge_category_csvs(args.years, categories=args.categories, output_path=combined_path)

    if args.plot:
        if args.output_plot:
            out = args.output_plot
        elif args.plot_type == "pareto":
            out = chart_path("sink_suitability_pareto.png")
        elif args.plot_type == "map":
            out = chart_path("sink_suitability_map.png")
        else:
            out = chart_path("sink_suitability_bivariate.png")

        if args.plot_type == "pareto":
            plot_pareto_panels(
                args.years, out,
                combined_path=combined_path,
                categories=args.categories,
            )
        elif args.plot_type == "bivariate":
            plot_bivariate_panels(
                args.years,
                out,
                combined_path=combined_path,
                categories=args.categories,
            )
        else:
            plot_map_panels(
                args.years,
                out,
                combined_path=combined_path,
                categories=args.categories,
            )

        summary_out = args.output_summary or data_path("sink_suitability_summary.csv")
        write_plot_summary_csv(
            args.years,
            summary_out,
            categories=args.categories,
            combined_path=combined_path,
        )


if __name__ == "__main__":
    main()
