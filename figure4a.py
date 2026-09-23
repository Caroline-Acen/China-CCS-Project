"""
Figure 4a: Seismic Screening Index (SSI) map
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

from config import NE_ADMIN0, NE_ADMIN1, chart_path, manuscript_font_bundle
from seismic_risk import SSI_COLORBAR_LABEL, srf_listed_colormap, style_srf_colorbar


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
        csv_path: Path to CSV with lat, lon, final_risk_score, storage_potential columns
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
    required = {"lat", "lon", "final_risk_score", "storage_potential"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    df = df[list(required)].copy()
    df = df.dropna(subset=["lat", "lon", "final_risk_score", "storage_potential"])
    df["final_risk_score"] = pd.to_numeric(df["final_risk_score"], errors="coerce")
    df["storage_potential"] = pd.to_numeric(df["storage_potential"], errors="coerce")
    df = df.dropna(subset=["final_risk_score", "storage_potential"]).drop_duplicates()

    # Compute storage potential totals by risk class (convert Mt -> Gt)
    low_total_gt = df.loc[df["final_risk_score"] <= 3, "storage_potential"].sum() / 1000.0
    moderate_total_gt = df.loc[
        (df["final_risk_score"] > 3) & (df["final_risk_score"] <= 6),
        "storage_potential",
    ].sum() / 1000.0
    high_total_gt = df.loc[
        (df["final_risk_score"] > 6) & (df["final_risk_score"] <= 9),
        "storage_potential",
    ].sum() / 1000.0
    very_high_total_gt = df.loc[df["final_risk_score"] > 9, "storage_potential"].sum() / 1000.0

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

    # Low (0–3) greens → Moderate (3–6) blues → High / Very High (>6) reds
    cmap, norm, boundaries = srf_listed_colormap(vmin, vmax)

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
    mf = manuscript_font_bundle(fig.get_figwidth())

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

    # Reserve right/bottom margins: colorbar label gap + room for risk summary
    fig.subplots_adjust(left=0.02, right=0.80, top=0.98, bottom=0.24)
    cax = fig.add_axes([0.83, 0.28, 0.018, 0.55])

    cbar = plt.colorbar(
        sm,
        cax=cax,
        orientation="vertical",
        boundaries=boundaries,
        ticks=ticks,
        spacing="proportional",
        drawedges=True,
    )
    style_srf_colorbar(
        fig,
        cbar,
        SSI_COLORBAR_LABEL,
        mf["colorbar_label"],
        round(mf["colorbar_tick"] * 0.95, 1),
        label_dx=0.04,
    )

    # Risk summary kept clear of the map (two lines)
    legend_text = (
        f"0≤SSI≤3 (Low Risk: {low_total_gt:,.2f} Gt) | "
        f"3<SSI≤6 (Moderate Risk: {moderate_total_gt:,.2f} Gt) |\n"
        f"6<SSI≤9 (High Risk: {high_total_gt:,.2f} Gt) | "
        f"SSI>9 (Very High Risk: {very_high_total_gt:,.2f} Gt)"
    )
    fig.text(0.5, 0.055, legend_text, ha="center", va="center", fontsize=mf["legend"])

    ax.axis("off")
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print(f"[OK] Map saved to {output_path}")


# Run function
plot_china_score_traces(
    csv_path="./data/Risk_Assessment/final_seismic_risk_factor_base_case.csv",
    admin0_path=str(NE_ADMIN0),
    admin1_path=str(NE_ADMIN1),
    output_path=chart_path("seismic_risk_factor_map.png"),
)
