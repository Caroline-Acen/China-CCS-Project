"""
Figure 3a: Seismic Risk Factor Map
China map showing locally-connected risk traces from seismic risk CSV data.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from sklearn.neighbors import NearestNeighbors


def plot_china_score_traces(
    csv_path,
    admin0_path,
    admin1_path,
    output_path="china_score_traces.png",
    target_crs="ESRI:102012",
    link_dist_m=60_000,
    k_neighbors=4,
    line_width=0.18,
    line_alpha=1.0,
    point_size=1.2,
    point_alpha=1.0,
    dpi=600,
):
    """
    Plot China map with locally-connected risk traces from point CSV.

    Parameters:
        csv_path: Path to CSV with lat, lon, final_risk_score columns
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
        target_crs: Target coordinate reference system
        link_dist_m: Maximum distance (meters) for connecting nearby points
        k_neighbors: Number of nearest neighbors to consider
        line_width: Width of trace lines
        point_size: Size of scatter points
        dpi: Output resolution
    """

    # Load and validate CSV data
    df = pd.read_csv(csv_path, low_memory=False)
    required = {"lat", "lon", "final_risk_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    df = df[list(required)].copy()
    df = df.dropna(subset=["lat", "lon", "final_risk_score"])
    df["final_risk_score"] = pd.to_numeric(df["final_risk_score"], errors="coerce")
    df = df.dropna(subset=["final_risk_score"]).drop_duplicates()

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"]),
        crs="EPSG:4326",
    ).to_crs(target_crs)

    # Determine color scale range
    vmin = 0
    vmax = int(np.ceil(float(gdf["final_risk_score"].max())))

    # Load and filter boundaries to China only
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

    # Clip provinces to China outline
    try:
        provinces = gpd.overlay(provinces, china0[["geometry"]], how="intersection")
    except Exception:
        provinces = provinces.copy()
        provinces["geometry"] = provinces.geometry.buffer(0)
        china_fixed = china0.copy()
        china_fixed["geometry"] = china_fixed.geometry.buffer(0)
        provinces = gpd.overlay(provinces, china_fixed[["geometry"]], how="intersection")

    # Setup discrete colormap
    boundaries = np.arange(vmin, vmax + 2)
    cmap = plt.get_cmap("turbo", len(boundaries) - 1)
    norm = BoundaryNorm(boundaries, cmap.N)

    # Build neighbor-based line segments
    xy = np.column_stack([gdf.geometry.x.values, gdf.geometry.y.values])
    scores = gdf["final_risk_score"].values
    segments = []
    seg_vals = []

    n = len(xy)
    if n >= 2:
        k = min(k_neighbors, n)
        nn = NearestNeighbors(n_neighbors=k, algorithm="auto")
        nn.fit(xy)
        dists, idxs = nn.kneighbors(xy)

        for i in range(n):
            for jpos in range(1, k):
                j = idxs[i, jpos]
                if dists[i, jpos] <= link_dist_m:
                    segments.append([xy[i], xy[j]])
                    seg_vals.append(scores[i])

    # Create figure
    fig, ax = plt.subplots(figsize=(11, 7))

    # Plot boundaries
    china0.boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
    provinces.boundary.plot(ax=ax, color="0.75", linewidth=0.6, zorder=2)

    # Set map extent
    minx, miny, maxx, maxy = china0.total_bounds
    pad_x = (maxx - minx) * 0.02
    pad_y = (maxy - miny) * 0.02
    ax.set_xlim(minx - pad_x, maxx + pad_x)
    ax.set_ylim(miny - pad_y, maxy + pad_y)

    # Plot line segments
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

    # Plot scatter points
    ax.scatter(
        gdf.geometry.x.values,
        gdf.geometry.y.values,
        c=scores,
        cmap=cmap,
        norm=norm,
        s=point_size,
        linewidths=0,
        alpha=point_alpha,
        zorder=4,
        rasterized=True,
    )

    # Add colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    # Determine tick spacing
    if vmax <= 12:
        step = 1
    elif vmax <= 25:
        step = 2
    elif vmax <= 50:
        step = 5
    else:
        step = 10

    ticks = list(np.arange(vmin, vmax + 1, step))
    if ticks[0] != vmin:
        ticks = [vmin] + ticks
    if ticks[-1] != vmax:
        ticks.append(vmax)

    # Auto-scale colorbar width
    n_bins = len(boundaries) - 1
    width = min(0.92, max(0.32, 0.04 * n_bins))
    left = 0.5 - width / 2
    cax = ax.inset_axes([left, -0.085, width, 0.028])

    cbar = plt.colorbar(
        sm,
        cax=cax,
        orientation="horizontal",
        boundaries=boundaries,
        ticks=ticks,
        spacing="proportional",
        drawedges=True,
    )
    cbar.set_label("Seismic risk factor (SRF)")
    cbar.outline.set_linewidth(0.6)
    try:
        cbar.dividers.set_linewidth(0.4)
    except Exception:
        pass

    ax.axis("off")
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Map saved to {output_path}")


# Run function
plot_china_score_traces(
    csv_path="./data/Risk_Assessment/final_seismic_risk_factor_base_case.csv",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="seismic_risk_factor_map.png"
)
