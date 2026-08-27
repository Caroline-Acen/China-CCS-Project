"""
Figure 20 Supplementary: Storage Potential and Industries Map
Maps showing storage potential locations overlaid with different industry categories.
"""

from pathlib import Path

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from shapely.ops import unary_union

from config import DSA_PIPELINE, EOR_PIPELINE, INDUSTRIES_CSV, NE_ADMIN0, NE_ADMIN1, CHARTS_DIR, manuscript_font_bundle


def plot_storage_potential_and_industries(
    dsa_pipeline, 
    eor_pipeline,
    industries_csv, 
    admin0_path, 
    admin1_path,
    output_prefix="storage_industry_plot"
):
    """
    Plot storage potential and industry locations with 250km buffer overlap analysis.

    Parameters:
        dsa_pipeline: Path to DSA Pipeline Excel file
        eor_pipeline: Path to EOR Pipeline Excel file
        industries_csv: Path to industries CSV
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_prefix: Prefix for output files
        
    Creates 4 plots with different industry category pairs.
    """

    print("Loading data...")
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

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

    # Load Industries in chunks (full file is ~97M rows)
    chunk_list = []
    industry_keys = list(categories.keys())
    for chunk in pd.read_csv(industries_csv, chunksize=250_000, low_memory=False):
        filtered = chunk[chunk[industry_keys].sum(axis=1) > 0]
        if not filtered.empty:
            chunk_list.append(filtered)
    ind_df = pd.concat(chunk_list, ignore_index=True)

    # Filter rows where all categories are zero
    ind_df = ind_df[ind_df[list(categories.keys())].sum(axis=1) > 0]
    print(f"Filtered industries: {len(ind_df)} rows kept")

    # Coordinate handling
    x_min, x_max = ind_df["X"].min(), ind_df["X"].max()
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

    # Shared color scale
    max_val = base_gdf[list(categories.keys())].max().max()
    bins = [0, 1, 2, 3, 4, 5, max_val]
    labels = ["1", "2", "3", "4", "5", f">{5}"]
    cmap = plt.get_cmap("YlOrRd", len(labels))
    norm = mcolors.BoundaryNorm(boundaries=bins, ncolors=len(labels), clip=False)

    # Legend handles
    legend_handles = []
    for cat, marker in categories.items():
        clean_label = cat.replace("19", "")
        legend_handles.append(Line2D([0], [0], marker=marker, color="w",
                                     markerfacecolor="gray", markersize=12, label=clean_label))
    legend_handles.append(Line2D([0], [0], marker="o", color="w",
                                 markerfacecolor="green", markersize=12, label="DSA storage potential"))
    legend_handles.append(Line2D([0], [0], marker="o", color="w",
                                 markerfacecolor="#2166ac", markersize=12, label="EOR storage potential"))
    legend_handles.append(Patch(facecolor="none", edgecolor="#6a3d9a", hatch="///",
                                label="Industries and storage\npotential < 250km"))

    # Create plots for each category pair
    category_pairs = [
        (["TC19", "ME19"], "Transport and Manufacturing"),
        (["WF19", "AF19"], "Waste and Agriculture"), 
        (["OM19", "MC19"], "Other Industries and Mining"),
        (["PM19", "PP19"], "Paper and Power")
    ]

    for plot_idx, (cat_pair, plot_title) in enumerate(category_pairs):
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
        mf = manuscript_font_bundle(fig.get_figwidth())
        
        for ax_idx, cat in enumerate(cat_pair):
            ax = axes[ax_idx]
            
            # Base maps
            china.boundary.plot(ax=ax, color="black", linewidth=1)
            provinces.boundary.plot(ax=ax, color="gray", linewidth=0.5, alpha=0.7)

            # Industries for this category
            subset = base_gdf[base_gdf[cat] > 0]
            if not subset.empty:
                subset.plot(
                    ax=ax,
                    markersize=40,
                    marker=categories[cat],
                    c=subset[cat],
                    cmap=cmap,
                    norm=norm,
                    alpha=0.7
                )

                # Calculate overlap
                buffer_250km = subset.buffer(250_000)
                buffer_union = unary_union(buffer_250km)

                overlap_polygons = []
                for pt in pd.concat([dsa_gdf, eor_gdf]).geometry:
                    if buffer_union.contains(pt):
                        overlap_polygons.append(buffer_union)
                
                if overlap_polygons:
                    overlap_gdf = gpd.GeoDataFrame(
                        geometry=[unary_union(overlap_polygons)],
                        crs=target_crs
                    )
                    overlap_gdf.plot(ax=ax, facecolor="none", edgecolor="#6a3d9a",
                                     hatch="///", linewidth=0, alpha=0.5)

            # Storage potentials
            dsa_gdf.plot(ax=ax, color="green", markersize=25, alpha=0.7)
            eor_gdf.plot(ax=ax, color="#2166ac", markersize=25, alpha=0.8)

            ax.set_title(cat.replace("19", ""), fontsize=mf["title"], fontweight="bold")
            ax.axis("off")

        # Add legend and colorbar to last plot
        if plot_idx == len(category_pairs) - 1:
            axes[-1].legend(handles=legend_handles, fontsize=mf["legend"],
                           loc="lower left", bbox_to_anchor=(-0.2, 0.05), frameon=False)

            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=axes.ravel().tolist(),
                               orientation="vertical", fraction=0.025, pad=-0.18)
            cbar.set_ticks([(bins[i] + bins[i+1]) / 2 for i in range(len(labels))])
            cbar.set_ticklabels(labels)
            cbar.set_label("Number of Industries", fontsize=mf["colorbar_label"], fontweight="bold")
            cbar.ax.tick_params(labelsize=mf["colorbar_tick"])

        plt.tight_layout()
        output_path = output_prefix.parent / f"{output_prefix.name}_{plot_idx + 1}.png"
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"[OK] Plot {plot_idx + 1} saved to {output_path}")
        plt.close()

    print("All plots generated successfully!")


# Run function
plot_storage_potential_and_industries(
    dsa_pipeline=str(DSA_PIPELINE),
    eor_pipeline=str(EOR_PIPELINE),
    industries_csv=str(INDUSTRIES_CSV),
    admin0_path=str(NE_ADMIN0),
    admin1_path=str(NE_ADMIN1),
    output_prefix=str(CHARTS_DIR / "storage_industry"),
)
