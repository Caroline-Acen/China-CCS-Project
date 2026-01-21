"""
Figure 1 Supplementary: Storage Potential and Industries Map (Original)
China map showing all industry categories overlaid with storage potential locations.
Similar to figure4a but used for supplementary materials.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from shapely.ops import unary_union
from textwrap import wrap


def plot_storage_potential_and_industries(
    dsa_pipeline, 
    eor_pipeline,
    industries_csv, 
    admin0_path, 
    admin1_path,
    output_path="storage_industry_map_sup.png"
):
    """
    Plot storage potential and industry locations with 250km buffer overlap.

    Parameters:
        dsa_pipeline: Path to DSA Pipeline Excel file
        eor_pipeline: Path to EOR Pipeline Excel file
        industries_csv: Path to industries CSV
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
    """

    print("Loading data...")

    # Load Storage Potential
    dsa_df = pd.read_excel(dsa_pipeline)
    eor_df = pd.read_excel(eor_pipeline)

    dsa_gdf = gpd.GeoDataFrame(
        dsa_df,
        geometry=gpd.points_from_xy(dsa_df["X"], dsa_df["Y"]),
        crs="EPSG:4326"
    )
    eor_gdf = gpd.GeoDataFrame(
        eor_df,
        geometry=gpd.points_from_xy(eor_df["X"], eor_df["Y"]),
        crs="EPSG:4326"
    )
    
    # Load Industries
    ind_df = pd.read_csv(industries_csv, low_memory=False)

    # Industry category markers
    categories = {
        "TC19": "o",   # Transport
        "ME19": "s",   # Manufacturing
        "WF19": "^",   # Waste
        "AF19": "D",   # Agriculture
        "OM19": "P",   # Other industries
        "MC19": "X",   # Mining
        "PM19": "v",   # Paper
        "PP19": "*",   # Power
    }

    # Filter rows where all categories are zero
    ind_df = ind_df[ind_df[list(categories.keys())].sum(axis=1) > 0]
    print(f"Filtered industries: {len(ind_df)} rows kept")

    # Coordinate handling
    x_min = ind_df["X"].min()
    target_crs = "ESRI:102012"
    already_projected = (abs(x_min) > 1000)

    if already_projected:
        base_gdf = gpd.GeoDataFrame(
            ind_df,
            geometry=gpd.points_from_xy(ind_df["X"], ind_df["Y"]),
            crs=target_crs
        )
    else:
        base_gdf = gpd.GeoDataFrame(
            ind_df,
            geometry=gpd.points_from_xy(ind_df["X"], ind_df["Y"]),
            crs="EPSG:4326"
        )

    # Load shapefiles
    admin0 = gpd.read_file(admin0_path)
    admin1 = gpd.read_file(admin1_path)
    
    china = admin0[admin0["ADMIN"] == "China"]
    if "admin" in admin1.columns:
        provinces = admin1[admin1["admin"] == "China"]
    elif "ADM0NAME" in admin1.columns:
        provinces = admin1[admin1["ADM0NAME"] == "China"]
    else:
        provinces = gpd.GeoDataFrame(columns=admin1.columns, geometry=admin1.geometry.name)
    
    # Reproject all data
    china = china.to_crs(target_crs)
    provinces = provinces.to_crs(target_crs)
    dsa_gdf = dsa_gdf.to_crs(target_crs)
    eor_gdf = eor_gdf.to_crs(target_crs)
    base_gdf = base_gdf.to_crs(target_crs)

    # Compute 250km buffers
    print("Computing 250km buffers...")
    buffer_250km = base_gdf.buffer(250_000)
    buffer_union = unary_union(buffer_250km)

    overlap_polygons = []
    for pt in pd.concat([dsa_gdf, eor_gdf]).geometry:
        if buffer_union.contains(pt):
            overlap_polygons.append(buffer_union)
    
    overlap_gdf = gpd.GeoDataFrame(
        geometry=[unary_union(overlap_polygons)],
        crs=target_crs
    ) if overlap_polygons else gpd.GeoDataFrame(geometry=[], crs=target_crs)

    # Create figure
    print("Creating plot...")
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    china.boundary.plot(ax=ax, color="black", linewidth=1)
    provinces.boundary.plot(ax=ax, color="gray", linewidth=0.5, alpha=0.7)

    # Color scale
    max_val = base_gdf[list(categories.keys())].max().max()
    bins = [0, 1, 2, 3, 4, 5, max_val]
    labels = ["1", "2", "3", "4", "5", f">{5}"]
    cmap = plt.get_cmap("YlOrRd", len(labels))
    norm = mcolors.BoundaryNorm(boundaries=bins, ncolors=len(labels), clip=False)

    # Plot industries by category
    for cat, marker in categories.items():
        subset = base_gdf[base_gdf[cat] > 0]
        if subset.empty:
            continue
        subset.plot(
            ax=ax,
            markersize=40,
            marker=marker,
            c=subset[cat],
            cmap=cmap,
            norm=norm,
            alpha=0.7
        )

    # Plot storage potentials
    dsa_gdf.plot(ax=ax, color="green", markersize=25, alpha=0.7)
    eor_gdf.plot(ax=ax, color="gray", markersize=25, alpha=0.7)

    # Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.03, pad=0.02)
    cbar.set_ticks([(bins[i] + bins[i+1]) / 2 for i in range(len(labels))])
    cbar.set_ticklabels(labels)
    cbar.set_label("Number of Industries", fontsize=10)

    # Overlap region
    if not overlap_gdf.empty:
        overlap_gdf.plot(ax=ax, facecolor="none", edgecolor="blue",
                         hatch="///", linewidth=0, alpha=0.5)

    # Legend
    legend_elements = []
    for cat, marker in categories.items():
        legend_elements.append(Line2D([0], [0], marker=marker, color="w",
                                      markerfacecolor="orange", markersize=8,
                                      label=cat.replace("19", "")))
    legend_elements.extend([
        Patch(facecolor="green", edgecolor="black", label="DSA storage potential"),
        Patch(facecolor="gray", edgecolor="black", label="EOR storage potential"),
        Patch(facecolor="none", edgecolor="blue", hatch="///",
              label="Industries and storage potential < 250km"),
    ])
    labels_wrapped = ["\n".join(wrap(elem.get_label(), 30)) for elem in legend_elements]
    ax.legend(handles=legend_elements, labels=labels_wrapped,
              loc="lower left", fontsize=9, frameon=True)

    ax.axis("off")
    
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Map saved to {output_path}")


# Run function
plot_storage_potential_and_industries(
    dsa_pipeline="./data/DSA-Pipeline.xlsx",
    eor_pipeline="./data/EOR_Pipeline.xlsx",
    industries_csv="./data/industries.csv",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="storage_industry_map_sup.png"
)
