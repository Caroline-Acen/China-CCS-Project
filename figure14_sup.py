"""
Figure 14 Supplementary: Seismic Risk Score (SRS) Map
China map showing seismic risk scores from CAFD fault proximity data.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import ScalarMappable
from sklearn.neighbors import NearestNeighbors


def plot_china_srf_with_provinces(
    srf_csv,
    admin0_path,
    admin1_path,
    output_path="seismic_risk_china.png",
    target_crs="ESRI:102012",
    link_dist_m=60_000,
    line_width=0.18,
    line_alpha=1.0,
    point_size=1.2,
    point_alpha=1.0,
    dpi=600,
):
    """
    Plot China map with Seismic Risk Score (SRS) traces from fault proximity data.

    Parameters:
        srf_csv: Path to CSV with X, Y, SRF columns
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
        target_crs: Target coordinate reference system
        link_dist_m: Maximum distance for connecting nearby points
        dpi: Output resolution
    """

    # Load SRF data
    df = pd.read_csv(srf_csv, low_memory=False)
    required = {"X", "Y", "SRF"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    df = df[["X", "Y", "SRF"]].copy()
    df = df.dropna(subset=["X", "Y", "SRF"])
    df["SRF"] = pd.to_numeric(df["SRF"], errors="coerce")
    df = df.dropna(subset=["SRF"]).drop_duplicates()

    srf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["X"], df["Y"]),
        crs="EPSG:4326",
    ).to_crs(target_crs)

    vmin = 0
    vmax = int(np.ceil(float(srf["SRF"].max())))

    # Load and filter boundaries
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

    # Setup discrete colormap
    boundaries = np.arange(vmin, vmax + 2)
    cmap = plt.get_cmap("turbo", len(boundaries) - 1)
    norm = BoundaryNorm(boundaries, cmap.N)

    # Build line segments by SRF class
    segments = []
    seg_vals = []

    for srf_val, grp in srf.groupby("SRF"):
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
                    segments.append([xy[i], xy[j]])
                    seg_vals.append(srf_val)

    # Create figure
    fig, ax = plt.subplots(figsize=(11, 7))

    # Plot boundaries
    china0.boundary.plot(ax=ax, color="black", linewidth=1.0, zorder=1)
    provinces.boundary.plot(ax=ax, color="0.75", linewidth=0.6, zorder=2)

    # Set extent
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

    # Plot points
    ax.scatter(
        srf.geometry.x.values,
        srf.geometry.y.values,
        c=srf["SRF"].values,
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

    ticks = list(range(vmin, vmax + 1))
    n_bins = len(boundaries) - 1
    width = min(0.92, max(0.35, 0.04 * n_bins))
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
    cbar.set_label("Seismic risk score (SRS)")
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
plot_china_srf_with_provinces(
    srf_csv="./data/Risk_Assessment/CAFD400-SRF-5km-2.csv",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="seismic_risk_srs_map.png"
)
