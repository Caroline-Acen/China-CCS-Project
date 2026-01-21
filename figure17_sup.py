"""
Figure 17 Supplementary: Sensitivity Analysis - PGA Thresholds
Multi-panel maps showing seismic risk factor at different PGA exceedance thresholds.
"""

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from sklearn.neighbors import NearestNeighbors


def plot_china_score_traces_multiple(
    admin0_path,
    admin1_path,
    output_path="china_score_traces_multi.png",
    target_crs="ESRI:102012",
    link_dist_m=60_000,
    k_neighbors=4,
    line_width=0.18,
    line_alpha=1.0,
    point_size=1.2,
    point_alpha=1.0,
    dpi=600,
    ncols=3,
    **csv_paths,
):
    """
    Multi-panel China risk-trace maps (one subplot per CSV).

    Parameters:
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
        csv_paths: CSV paths as csv_path1=..., csv_path2=..., etc.
                   Each CSV must have lat, lon, final_risk_score columns
    """

    # Extract CSV paths in order
    def _extract_index(k):
        digits = "".join([c for c in k if c.isdigit()])
        return int(digits) if digits else 10**9

    csv_items = [(k, v) for k, v in csv_paths.items() if k.lower().startswith("csv_path")]
    if not csv_items:
        raise ValueError("No csv_path* arguments provided.")

    csv_items = sorted(csv_items, key=lambda kv: _extract_index(kv[0]))
    csv_files = [v for _, v in csv_items]
    n_panels = len(csv_files)

    ncols = max(1, int(ncols))
    nrows = int(np.ceil(n_panels / ncols))

    # Load boundaries
    admin0 = gpd.read_file(admin0_path).to_crs(target_crs)
    admin1 = gpd.read_file(admin1_path).to_crs(target_crs)

    # Find China
    china0 = None
    for col in ["ADMIN", "NAME", "SOVEREIGNT", "ADM0NAME", "NAME_EN", "COUNTRY"]:
        if col in admin0.columns:
            china0 = admin0[admin0[col].astype(str).str.strip().str.lower().eq("china")]
            if not china0.empty:
                break
    if china0 is None or china0.empty:
        raise ValueError("China not found in admin0 shapefile.")

    # Find provinces
    provinces = None
    for col in ["admin", "ADM0NAME", "ADMIN", "SOVEREIGNT", "COUNTRY"]:
        if col in admin1.columns:
            provinces = admin1[admin1[col].astype(str).str.strip().str.lower().eq("china")]
            if provinces is not None and not provinces.empty:
                break
    if provinces is None or provinces.empty:
        provinces = admin1.copy()

    # Clip provinces to China
    try:
        provinces = gpd.overlay(provinces, china0[["geometry"]], how="intersection")
    except Exception:
        provinces = provinces.copy()
        provinces["geometry"] = provinces.geometry.buffer(0)
        china_fixed = china0.copy()
        china_fixed["geometry"] = china_fixed.geometry.buffer(0)
        provinces = gpd.overlay(provinces, china_fixed[["geometry"]], how="intersection")

    # Set extent
    minx, miny, maxx, maxy = china0.total_bounds
    pad_x = (maxx - minx) * 0.02
    pad_y = (maxy - miny) * 0.02
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y)

    # Load all CSVs and compute global max
    gdfs = []
    global_max = 0.0

    for path in csv_files:
        df = pd.read_csv(path, low_memory=False)
        required = {"lat", "lon", "final_risk_score"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{path} missing columns {missing}")

        df = df[["lat", "lon", "final_risk_score"]].copy()
        df = df.dropna(subset=["lat", "lon", "final_risk_score"])
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

    # Setup colormap
    boundaries = np.arange(vmin, vmax + 2)
    cmap = plt.get_cmap("turbo", len(boundaries) - 1)
    norm = BoundaryNorm(boundaries, cmap.N)

    # Calculate ticks
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

    # Create figure
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(11 * ncols, 7 * nrows), constrained_layout=False)
    axes = np.array(axes).reshape(-1)
    fig.subplots_adjust(left=0.02, right=0.7, top=0.98, bottom=0.12, wspace=0.02, hspace=0.05)

    letters = "abcdefghijklmnopqrstuvwxyz"

    # Draw each panel
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
        ax.text(0.5, 0.995, lab, transform=ax.transAxes, ha="center", va="top", fontsize=14, fontweight="bold", zorder=10)

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
            lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=line_width, alpha=line_alpha, zorder=3, rasterized=True)
            lc.set_array(np.asarray(seg_vals))
            ax.add_collection(lc)

        ax.scatter(xy[:, 0], xy[:, 1], c=scores, cmap=cmap, norm=norm, s=point_size, linewidths=0, alpha=point_alpha, zorder=4, rasterized=True)

    # Shared colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    n_bins = len(boundaries) - 1
    width = min(0.92, max(0.32, 0.04 * n_bins))
    left = 0.5 - width / 2 - 0.1

    cax = fig.add_axes([left, 0.05, width, 0.018])
    cbar = plt.colorbar(sm, cax=cax, orientation="horizontal", boundaries=boundaries, ticks=ticks, spacing="proportional", drawedges=True)
    cbar.set_label("Seismic risk factor (SRF)")
    cbar.outline.set_linewidth(0.6)

    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved multi-panel figure to {output_path}")


# Run function
plot_china_score_traces_multiple(
    csv_path1="./data/Risk_Assessment/sensitivity_pga_0.5_risk.csv",
    csv_path2="./data/Risk_Assessment/sensitivity_pga_1_risk.csv",
    csv_path3="./data/Risk_Assessment/sensitivity_pga_2_risk.csv",
    csv_path4="./data/Risk_Assessment/sensitivity_pga_5_risk.csv",
    csv_path5="./data/Risk_Assessment/sensitivity_pga_63_risk.csv",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="sensitivity_pga.png"
)
