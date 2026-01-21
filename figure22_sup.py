"""
Figure 22 Supplementary: VS30 Map
China map showing VS30 (shear wave velocity at 30m depth) values.
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from sklearn.neighbors import NearestNeighbors
from matplotlib.collections import LineCollection


def plot_china_vs30_map(
    csv_path,
    admin0_path,
    admin1_path,
    output_path="vs30_map.png",
    target_crs="ESRI:102012",
    link_dist_m=60_000,
    k_neighbors=4,
    line_width=0.18,
    line_alpha=1.0,
    point_size=3.0,
    point_alpha=1.0,
    dpi=600,
):
    """
    Single-panel China VS30 map.

    Parameters:
        csv_path: Path to CSV with Long, Lat, Vs30 columns
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
        target_crs: Target coordinate reference system
        link_dist_m: Maximum distance for connecting nearby points
        point_size: Size of scatter points
        dpi: Output resolution
    """

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

    # Merge China + Taiwan geometries
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

    # Load VS30 data
    df = pd.read_csv(csv_path, low_memory=False)
    required = {"Long", "Lat", "Vs30"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns {missing}")

    df = df[["Long", "Lat", "Vs30"]].copy()
    df = df.dropna(subset=["Long", "Lat", "Vs30"])
    df["Vs30"] = pd.to_numeric(df["Vs30"], errors="coerce")
    df = df.dropna(subset=["Vs30"]).drop_duplicates()

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["Long"], df["Lat"]),
        crs="EPSG:4326",
    ).to_crs(target_crs)

    # Clip to China
    gdf = gdf[gdf.within(china_geom)].copy()

    # Setup colormap
    vmin = 0
    vmax = int(np.ceil(gdf["Vs30"].max() / 100) * 100)

    step = 100
    boundaries = np.arange(vmin, vmax + step, step)
    cmap = plt.get_cmap("viridis", len(boundaries) - 1)
    norm = BoundaryNorm(boundaries, cmap.N)

    # Ticks
    if vmax <= 500:
        tick_step = 100
    elif vmax <= 1000:
        tick_step = 200
    else:
        tick_step = 250
    ticks = list(np.arange(vmin, vmax + 1, tick_step))
    if ticks[-1] != vmax:
        ticks.append(vmax)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10), constrained_layout=False)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.12)

    # Plot boundaries - show China and Taiwan separately for distinct outline
    china0.boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
    if taiwan0 is not None and not taiwan0.empty:
        taiwan0.to_crs(target_crs).boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
    provinces.boundary.plot(ax=ax, color="0.75", linewidth=0.6, zorder=2)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")

    # Build segments
    xy = np.column_stack([gdf.geometry.x.values, gdf.geometry.y.values])
    scores = gdf["Vs30"].values

    segments = []
    seg_vals = []

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
        c=scores,
        cmap=cmap,
        norm=norm,
        s=point_size,
        linewidths=0,
        alpha=point_alpha,
        zorder=4,
        rasterized=True,
    )

    # Colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    n_bins = len(boundaries) - 1
    width = min(0.7, max(0.4, 0.05 * n_bins))
    left = 0.5 - width / 2

    cax = fig.add_axes([left, 0.05, width, 0.02])
    cbar = plt.colorbar(
        sm,
        cax=cax,
        orientation="horizontal",
        boundaries=boundaries,
        ticks=ticks,
        spacing="proportional",
        drawedges=True,
    )
    cbar.set_label("VS30 (m/s)", fontsize=12, fontweight="bold")
    cbar.outline.set_linewidth(0.6)
    try:
        cbar.dividers.set_linewidth(0.4)
    except Exception:
        pass

    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved VS30 map to {output_path}")


# Run function
plot_china_vs30_map(
    csv_path="./data/Risk_Assessment/VS30_China_5km_Interpolated.csv",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="vs30_map.png"
)
