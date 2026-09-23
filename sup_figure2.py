"""Supplementary Figure 2: Fault traces classified by age."""

import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from sklearn.neighbors import NearestNeighbors

from config import CAFD_FAULTS_CSV, NE_ADMIN0, NE_ADMIN1, chart_path, manuscript_font_bundle


def _resolve_first_existing_path(candidates):
    """Return the first existing path from a list of candidates."""
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return candidates[0] if candidates else None


def plot_china_faults_by_name(
    faults_csv,
    admin0_path,
    admin1_path,
    output_path="seismic_faults_by_name_map.png",
    target_crs="ESRI:102012",
    link_dist_m=60_000,
    line_width=0.22,
    line_alpha=0.92,
    point_size=1.0,
    point_alpha=0.55,
    target_column="AGE",
    top_n_categories=14,
    dpi=600,
):
    """
    Plot China map with fault traces classified by a selected fault attribute.

    Parameters:
        faults_csv: Path to CSV with X (lon), Y (lat), and optional class columns
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
        target_crs: Target coordinate reference system
        link_dist_m: Maximum distance for connecting nearby points
        target_column: Class column to color by (e.g., AGE or Fea_En)
        top_n_categories: If number of classes is large, keep top N and group others
        dpi: Output resolution
    """

    # Load fault data
    df = pd.read_csv(faults_csv, low_memory=False)

    required = {"X (lon)", "Y (lat)"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    # Resolve class field robustly (supports case differences like Fea_EN vs Fea_En).
    effective_column = target_column
    if effective_column not in df.columns:
        lower_lookup = {str(col).lower(): col for col in df.columns}
        if target_column.lower() in lower_lookup:
            effective_column = lower_lookup[target_column.lower()]
        else:
            # If requested class field does not exist, classify everything as Unknown.
            df[target_column] = "Unknown"
            effective_column = target_column

    df = df[["X (lon)", "Y (lat)", effective_column]].copy()
    df = df.dropna(subset=["X (lon)", "Y (lat)"])
    df["X (lon)"] = pd.to_numeric(df["X (lon)"], errors="coerce")
    df["Y (lat)"] = pd.to_numeric(df["Y (lat)"], errors="coerce")
    df[effective_column] = (
        df[effective_column]
        .fillna("Unknown")
        .astype(str)
        .str.strip()
        .replace("", "Unknown")
    )
    df = df.dropna(subset=["X (lon)", "Y (lat)"]).drop_duplicates()

    # Drop unknown classes completely as requested.
    df = df[df[effective_column].str.lower() != "unknown"].copy()
    if df.empty:
        raise ValueError(f"No valid rows left after dropping Unknown values in '{effective_column}'.")

    # Keep top classes if there are too many unique labels.
    counts = df[effective_column].value_counts()
    all_names = counts.index.tolist()
    if len(all_names) > top_n_categories:
        keep = set(all_names[:top_n_categories])
        df["FaultClass"] = np.where(df[effective_column].isin(keep), df[effective_column], "Other")
    else:
        df["FaultClass"] = df[effective_column]

    faults = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["X (lon)"], df["Y (lat)"]),
        crs="EPSG:4326",
    ).to_crs(target_crs)

    # Load boundaries
    if not admin0_path or not admin1_path:
        raise ValueError("Provide both admin0_path and admin1_path.")

    admin0 = gpd.read_file(admin0_path).to_crs(target_crs)
    admin1 = gpd.read_file(admin1_path).to_crs(target_crs)

    # Find China in admin0
    china0 = None
    for col in ["ADMIN", "NAME", "SOVEREIGNT", "ADM0NAME", "NAME_EN", "COUNTRY"]:
        if col in admin0.columns:
            china0 = admin0[admin0[col].astype(str).str.strip().str.lower().eq("china")]
            if not china0.empty:
                break
    if china0 is None or china0.empty:
        raise ValueError("China not found in admin0.")

    # Find China provinces
    provinces = None
    for col in ["admin", "ADM0NAME", "ADMIN", "SOVEREIGNT", "COUNTRY"]:
        if col in admin1.columns:
            provinces = admin1[admin1[col].astype(str).str.strip().str.lower().eq("china")]
            if provinces is not None and not provinces.empty:
                break
    if provinces is None or provinces.empty:
        provinces = admin1.copy()

    # Clip provinces to China outline
    try:
        provinces = gpd.overlay(provinces, china0[["geometry"]], how="intersection")
    except Exception:
        provinces = provinces.copy()
        provinces["geometry"] = provinces.geometry.buffer(0)
        china_fixed = china0.copy()
        china_fixed["geometry"] = china_fixed.geometry.buffer(0)
        provinces = gpd.overlay(provinces, china_fixed[["geometry"]], how="intersection")

    # Color map for categorical classes
    class_names = list(faults["FaultClass"].value_counts().index)
    n_classes = len(class_names)
    colors = plt.cm.tab20(np.linspace(0, 1, max(n_classes, 2)))
    class_to_color = {name: colors[i] for i, name in enumerate(class_names)}

    # Build line segments class by class
    class_segments = {name: [] for name in class_names}

    for fault_class, grp in faults.groupby("FaultClass"):
        xy = np.column_stack([grp.geometry.x.values, grp.geometry.y.values])
        n = len(xy)
        if n < 2:
            continue

        k = min(4, n)
        nn = NearestNeighbors(n_neighbors=k, algorithm="auto")
        nn.fit(xy)
        dists, idxs = nn.kneighbors(xy)

        for i in range(n):
            for jpos in range(1, k):
                j = idxs[i, jpos]
                if dists[i, jpos] <= link_dist_m:
                    class_segments[fault_class].append([xy[i], xy[j]])

    # Create figure
    fig, ax = plt.subplots(figsize=(11, 7.4))
    mf = manuscript_font_bundle(fig.get_figwidth())

    china0.boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
    provinces.boundary.plot(ax=ax, color="0.75", linewidth=0.6, zorder=2)

    # Set extent
    minx, miny, maxx, maxy = china0.total_bounds
    pad_x = (maxx - minx) * 0.02
    pad_y = (maxy - miny) * 0.02
    ax.set_xlim(minx - pad_x, maxx + pad_x)
    ax.set_ylim(miny - pad_y, maxy + pad_y)

    # Draw categorized line traces
    for cname, segs in class_segments.items():
        if not segs:
            continue
        lc = LineCollection(
            segs,
            colors=[class_to_color[cname]],
            linewidths=line_width,
            alpha=line_alpha,
            zorder=3,
            rasterized=True,
        )
        ax.add_collection(lc)

    # Draw points with the same class color
    for cname, grp in faults.groupby("FaultClass"):
        ax.scatter(
            grp.geometry.x.values,
            grp.geometry.y.values,
            c=[class_to_color[cname]],
            s=point_size,
            linewidths=0,
            alpha=point_alpha,
            zorder=4,
            rasterized=True,
        )

    # Reserve space for an external legend so it does not cover the map.
    fig.subplots_adjust(right=0.77)

    # Legend (single column for readability with long names)
    legend_order = faults["FaultClass"].value_counts().index.tolist()
    handles = [
        Line2D([0], [0], color=class_to_color[name], lw=2.0, label=name)
        for name in legend_order
    ]
    legend_title = f"{target_column} class"
    if effective_column.upper() == "AGE":
        legend_title = "Fault Ages"
    elif effective_column.lower() == "fea_en":
        legend_title = "Fault Features"

    ax.legend(
        handles=handles,
        title=legend_title,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
        fontsize=mf["legend"],
        title_fontsize=mf["title"],
        frameon=True,
        ncol=1,
    )

    ax.axis("off")
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Map saved to {output_path}")


if __name__ == "__main__":
    faults_csv_path = _resolve_first_existing_path(
        [
            str(CAFD_FAULTS_CSV),
            "./data/Risk_Assessment/CAFD400_V2023_1-Reprojected-5km-Intersection-Aggregate-GPSG-4326-Geom.csv",
        ]
    )
    plot_china_faults_by_name(
        faults_csv=faults_csv_path,
        admin0_path=str(NE_ADMIN0),
        admin1_path=str(NE_ADMIN1),
        target_column="AGE",
        output_path=chart_path("seismic_faults_by_name_map_v2.png"),
    )
