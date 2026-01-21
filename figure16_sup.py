"""
Figure 16 Supplementary: Multi-Panel PGA Maps
China maps showing Peak Ground Acceleration at different exceedance probabilities.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from sklearn.neighbors import NearestNeighbors


def plot_china_score_traces_multiple_pga(
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
    Multi-panel China PGA trace maps.

    Parameters:
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
        target_crs: Target coordinate reference system
        ncols: Number of columns in subplot grid
        csv_paths: CSV paths as csv_path1=..., csv_path2=..., etc.
                   Each CSV must have X, Y, 'PE50a0.5% (g)' columns

    Features:
        - China outline includes Taiwan
        - Points clipped to China boundary
        - Shared colorbar at bottom
        - Panels labeled (a), (b), (c)...
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
    for col in ["ADMIN", "NAME", "SOVEREIGNT", "ADM0NAME", "NAME_EN", "FORMAL_EN", "COUNTRY"]:
        if col in admin0.columns:
            china0 = admin0[admin0[col].astype(str).str.strip().str.lower().eq("china")]
            if not china0.empty:
                break
    if china0 is None or china0.empty:
        raise ValueError("China not found in admin0 shapefile.")

    # Find Taiwan
    taiwan0 = None
    patterns = ["taiwan", "taipei", "formosa", "chinese taipei"]
    for col in ["ADMIN", "NAME", "SOVEREIGNT", "ADM0NAME", "NAME_EN", "FORMAL_EN"]:
        if col in admin0.columns:
            s = admin0[col].astype(str).str.lower()
            for p in patterns:
                hit = admin0[s.str.contains(p, na=False)]
                if hit is not None and not hit.empty:
                    taiwan0 = hit
                    break
            if taiwan0 is not None and not taiwan0.empty:
                break

    # Merge China + Taiwan
    china_geom = china0.unary_union
    if taiwan0 is not None and not taiwan0.empty:
        china_geom = china_geom.union(taiwan0.unary_union)

    china_plot = gpd.GeoDataFrame(geometry=[china_geom], crs=target_crs)

    # Find provinces
    provinces = None
    for col in ["admin", "ADM0NAME", "ADMIN", "SOVEREIGNT", "COUNTRY"]:
        if col in admin1.columns:
            provinces = admin1[admin1[col].astype(str).str.strip().str.lower().eq("china")]
            if provinces is not None and not provinces.empty:
                break
    if provinces is None or provinces.empty:
        provinces = admin1.copy()

    # Clip provinces
    try:
        provinces = gpd.overlay(provinces, china_plot, how="intersection")
    except Exception:
        provinces = provinces.copy()
        provinces["geometry"] = provinces.geometry.buffer(0)
        china_fix = china_plot.copy()
        china_fix["geometry"] = china_fix.geometry.buffer(0)
        provinces = gpd.overlay(provinces, china_fix, how="intersection")

    # Set extent
    minx, miny, maxx, maxy = china_plot.total_bounds
    pad_x = (maxx - minx) * 0.02
    pad_y = (maxy - miny) * 0.02
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y)

    # Load all CSVs and compute global max PGA
    accel_col = "PE50a0.5% (g)"
    gdfs = []
    global_max = 0.0

    for path in csv_files:
        df = pd.read_csv(path, low_memory=False)

        required = {"X", "Y", accel_col}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{path} missing columns {missing}")

        df = df[["X", "Y", accel_col]].copy()
        df = df.dropna(subset=["X", "Y", accel_col])
        df[accel_col] = pd.to_numeric(df[accel_col], errors="coerce")
        df = df.dropna(subset=[accel_col]).drop_duplicates()

        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(df["X"], df["Y"]),
            crs="EPSG:4326",
        ).to_crs(target_crs)

        gdf = gdf[gdf.within(china_geom)].copy()
        gdfs.append(gdf)
        
        if len(gdf):
            global_max = max(global_max, float(gdf[accel_col].max()))

    vmin = 0.0
    vmax = float(global_max)

    # Set bin step for colormap
    if vmax <= 0.5:
        bin_step = 0.05
    elif vmax <= 1.5:
        bin_step = 0.1
    elif vmax <= 3.0:
        bin_step = 0.2
    else:
        bin_step = 0.5

    if vmax <= 0:
        boundaries = np.array([0.0, 1e-6])
    else:
        boundaries = np.arange(vmin, vmax + bin_step, bin_step)
        if boundaries[-1] < vmax:
            boundaries = np.append(boundaries, vmax)

    cmap = plt.get_cmap("turbo", max(2, len(boundaries) - 1))
    norm = BoundaryNorm(boundaries, cmap.N)

    # Calculate ticks
    max_ticks = 12
    tick_step = bin_step
    while (vmax / tick_step) > max_ticks:
        tick_step *= 2

    ticks = list(np.arange(0, vmax + 1e-12, tick_step))
    if not ticks:
        ticks = [0.0, vmax]
    elif ticks[-1] < vmax:
        ticks.append(vmax)

    # Create figure
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(11 * ncols, 7 * nrows),
        constrained_layout=False,
    )
    axes = np.array(axes).reshape(-1)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.12, wspace=0.02, hspace=0.02)

    letters = "abcdefghijklmnopqrstuvwxyz"

    # Draw each panel
    for i, ax in enumerate(axes):
        if i >= n_panels:
            ax.axis("off")
            continue

        gdf = gdfs[i]

        # Plot boundaries
        china_plot.boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
        provinces.boundary.plot(ax=ax, color="0.75", linewidth=0.6, zorder=2)

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.axis("off")

        # Panel label
        lab = f"({letters[i]})" if i < len(letters) else f"({i+1})"
        ax.text(
            0.5, 0.995, lab,
            transform=ax.transAxes,
            ha="center", va="top",
            fontsize=14, fontweight="bold",
            zorder=10,
        )

        if gdf.empty:
            continue

        xy = np.column_stack([gdf.geometry.x.values, gdf.geometry.y.values])
        pga = gdf[accel_col].values

        # Build segments
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
                        seg_vals.append(pga[p])

        # Plot segments
        if segments:
            lc = LineCollection(
                segments,
                cmap=cmap,
                norm=norm,
                linewidths=line_width,
                alpha=line_alpha,
                zorder=3,
                rasterized=True,
            )
            lc.set_array(np.asarray(seg_vals))
            ax.add_collection(lc)

        # Plot points
        ax.scatter(
            xy[:, 0], xy[:, 1],
            c=pga,
            cmap=cmap,
            norm=norm,
            s=point_size,
            linewidths=0,
            alpha=point_alpha,
            zorder=4,
            rasterized=True,
        )

    # Add shared colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    n_bins = len(boundaries) - 1
    width = min(0.92, max(0.32, 0.03 * n_bins))
    left = 0.5 - width / 2

    cax = fig.add_axes([left, 0.05, width, 0.018])
    cbar = plt.colorbar(
        sm,
        cax=cax,
        orientation="horizontal",
        boundaries=boundaries,
        ticks=ticks,
        spacing="proportional",
        drawedges=True,
    )
    cbar.set_label("Peak Ground Acceleration (g)")
    cbar.outline.set_linewidth(0.6)
    try:
        cbar.dividers.set_linewidth(0.4)
    except Exception:
        pass

    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved multi-panel figure to {output_path}")


# Run function
plot_china_score_traces_multiple_pga(
    csv_path1="./data/Risk_Assessment/PGA-PE50a0.5-IDW-ESRI-102012-5km-2.csv",
    csv_path2="./data/Risk_Assessment/PGA-PE50a1-IDW-ESRI-102012-5km.csv",
    csv_path3="./data/Risk_Assessment/PGA-PE50a2-IDW-ESRI-102012-5km.csv",
    csv_path4="./data/Risk_Assessment/PGA-PE50a5-IDW-ESRI-102012-5km.csv",
    csv_path5="./data/Risk_Assessment/PGA-PE50a10-IDW-ESRI-102012-5km.csv",
    csv_path6="./data/Risk_Assessment/PGA-PE50a63-IDW-ESRI-102012-5km.csv",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="pga_multi_panel.png"
)
