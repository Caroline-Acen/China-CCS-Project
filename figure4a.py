"""
Figure 4a: Storage Potential and Industries Combined Map
China map showing all industry categories overlaid with storage potential locations.
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
import os
import time

from config import (
    DSA_PIPELINE,
    EOR_PIPELINE,
    INDUSTRIES_CSV,
    NE_ADMIN0,
    NE_ADMIN1,
    chart_path,
    data_path,
    manuscript_font_bundle,
)


def plot_storage_potential_and_industries(
    dsa_pipeline, 
    eor_pipeline,
    industries_csv, 
    admin0_path, 
    admin1_path,
    output_path="storage_industry_map.png"
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

    t0 = time.perf_counter()

    def log_step(message):
        elapsed = time.perf_counter() - t0
        print(f"[{elapsed:8.2f}s] {message}")

    log_step("Loading data...")

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
    

    # Load Industries in chunks to avoid memory issues
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
    # Only keep rows where any category column is > 0
    chunk_list = []
    log_step("Reading industries CSV in chunks...")
    for idx, chunk in enumerate(pd.read_csv(industries_csv, low_memory=False, chunksize=100000), start=1):
        filtered = chunk[chunk[list(categories.keys())].sum(axis=1) > 0]
        chunk_list.append(filtered)
        if idx % 20 == 0:
            log_step(f"Processed {idx} chunks; retained rows so far: {sum(len(df) for df in chunk_list):,}")
    ind_df = pd.concat(chunk_list, ignore_index=True)

    # Filter rows where all categories are zero
    ind_df = ind_df[ind_df[list(categories.keys())].sum(axis=1) > 0]
    log_step(f"Filtered industries: {len(ind_df):,} rows kept")

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


    # Compute 250km proximity region around storage points
    log_step("Computing 250km buffers...")
    storage_points = gpd.GeoDataFrame(
        pd.concat([dsa_gdf, eor_gdf], ignore_index=True),
        geometry="geometry",
        crs=target_crs,
    )
    storage_union_buffer = unary_union(storage_points.geometry.buffer(250_000))
    log_step("Storage buffer union built; selecting industries within 250km...")

    # Chunked intersects avoids creating very large spatial-join intermediates
    chunk_size = 250_000
    mask_parts = []
    total_rows = len(base_gdf)
    for start in range(0, total_rows, chunk_size):
        end = min(start + chunk_size, total_rows)
        chunk_mask = base_gdf.geometry.iloc[start:end].intersects(storage_union_buffer)
        mask_parts.append(chunk_mask)
        processed = end
        if processed % 1_000_000 == 0 or processed == total_rows:
            log_step(f"Proximity scan progress: {processed:,}/{total_rows:,} rows")

    within_250km_mask = pd.concat(mask_parts).sort_index()
    industries_within_250km = base_gdf[within_250km_mask.values]
    unique_industries_within_250km = len(industries_within_250km)
    log_step(f"Unique industry rows within 250km: {unique_industries_within_250km:,}")

    # Calculate industry type counts and percentages
    industry_types = {
        "TC19": "TC",
        "ME19": "ME",
        "WF19": "WF",
        "AF19": "AF",
        "OM19": "OM",
        "MC19": "MC",
        "PM19": "PM",
        "PP19": "PP",
    }
    industry_keys = list(industry_types.keys())
    industry_counts = {k: (industries_within_250km[k] > 0).sum() for k in industry_keys}
    total_industry_types_within_250km = sum(industry_counts.values())
    industry_percentages = {
        k: (v / total_industry_types_within_250km * 100 if total_industry_types_within_250km > 0 else 0)
        for k, v in industry_counts.items()
    }

    # Use the same >0 category-membership logic for totals so category sums align with summary totals.
    total_category_assignments_within_250km = int((industries_within_250km[industry_keys] > 0).sum().sum())
    total_category_assignments_all_industries = int((base_gdf[industry_keys] > 0).sum().sum())
    category_assignments_within_250km_percentage = (
        total_category_assignments_within_250km / total_category_assignments_all_industries * 100
        if total_category_assignments_all_industries > 0 else 0
    )
    log_step(f"TOTAL_CATEGORY_ASSIGNMENTS_WITHIN_250KM: {total_category_assignments_within_250km:,}")

    multi_category_rows = int(((industries_within_250km[industry_keys] > 0).sum(axis=1) > 1).sum())
    log_step(f"Rows with multiple >0 categories within 250km: {multi_category_rows:,}")

    # Prepare DataFrame and save CSV (tabular output lives under data/)
    csv_path = data_path("storage_industry_map_industry_stats.csv")
    industry_stats_df = pd.DataFrame(
        {
            "Industry Type": [industry_types[k] for k in industry_types.keys()],
            "Count": [industry_counts[k] for k in industry_types.keys()],
            "Percentage": [round(industry_percentages[k], 2) for k in industry_types.keys()],
        }
    )
    summary_df = pd.DataFrame(
        {
            "Industry Type": [
                "TOTAL_CATEGORY_ASSIGNMENTS_WITHIN_250KM",
                "TOTAL_CATEGORY_ASSIGNMENTS_ALL_INDUSTRIES",
                "PERCENTAGE_CATEGORY_ASSIGNMENTS_WITHIN_250KM",
            ],
            "Count": [
                total_category_assignments_within_250km,
                total_category_assignments_all_industries,
                "",
            ],
            "Percentage": [
                round(category_assignments_within_250km_percentage, 2),
                100.0,
                round(category_assignments_within_250km_percentage, 2),
            ],
        }
    )
    industry_stats_df = pd.concat([industry_stats_df, summary_df], ignore_index=True)
    industry_stats_df.to_csv(csv_path, index=False)
    log_step(f"[OK] Industry stats CSV saved to {csv_path}")

    # Overlap polygons for plotting
    overlap_geom = storage_union_buffer
    overlap_gdf = gpd.GeoDataFrame(geometry=[overlap_geom], crs=target_crs)

    # Create figure
    log_step("Creating plot...")
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    mf = manuscript_font_bundle(fig.get_figwidth())

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
    cbar.set_label("Number of Industries", fontsize=mf["colorbar_label"])

    # Overlap region
    if not overlap_gdf.empty:
        overlap_gdf.plot(ax=ax, facecolor="none", edgecolor="blue",
                         hatch="///", linewidth=0, alpha=0.5)

    # Legend
    legend_elements = []
    for cat, marker in categories.items():
        legend_elements.append(Line2D([0], [0], marker=marker, color="w",
                                      markerfacecolor="orange", markersize=mf["legend_marker"],
                                      label=cat.replace("19", "")))
    legend_elements.extend([
        Patch(facecolor="green", edgecolor="black", label="DSA storage potential"),
        Patch(facecolor="gray", edgecolor="black", label="EOR storage potential"),
        Patch(facecolor="none", edgecolor="blue", hatch="///",
              label="Industries and storage potential < 250km"),
    ])
    labels_wrapped = ["\n".join(wrap(elem.get_label(), 30)) for elem in legend_elements]
    ax.legend(handles=legend_elements, labels=labels_wrapped,
              loc="lower left", fontsize=mf["legend"], frameon=True)

    ax.axis("off")
    
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    log_step(f"[OK] Map saved to {output_path}")


# Run function
plot_storage_potential_and_industries(
    dsa_pipeline=str(DSA_PIPELINE),
    eor_pipeline=str(EOR_PIPELINE),
    industries_csv=str(INDUSTRIES_CSV),
    admin0_path=str(NE_ADMIN0),
    admin1_path=str(NE_ADMIN1),
    output_path=chart_path("storage_industry_map.png"),
)
