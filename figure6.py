"""
Figure 6: Storage Potential and Gas Fields Map
China map showing storage potential locations (DSA/EOR) and gas field locations with 250km buffer overlap.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from textwrap import wrap
from shapely.ops import unary_union


def plot_storage_potential_and_gas_fields(
    dsa_pipeline, 
    eor_pipeline,
    china_oil_gas_csv, 
    admin0_path, 
    admin1_path,
    output_path="storage_gas_map.png"
):
    """
    Plot China map with storage potential and gas field locations.

    Parameters:
        dsa_pipeline: Path to DSA Pipeline Excel file
        eor_pipeline: Path to EOR Pipeline Excel file
        china_oil_gas_csv: Path to gas fields CSV
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path

    Markers:
        - Green: DSA storage potential
        - Orange: EOR storage potential
        - Red: Gas fields
        - Blue hatched: Areas where storage and gas fields are within 250km
    """
    
    # Load Excel files (Storage Potential)
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
    
    # Load Gas Fields CSV
    gas_df = pd.read_csv(china_oil_gas_csv, low_memory=False)
    gas_gdf = gpd.GeoDataFrame(
        gas_df,
        geometry=gpd.points_from_xy(gas_df["X"], gas_df["Y"]),
        crs="EPSG:4326"
    )
    
    # Load Natural Earth shapefiles
    admin0 = gpd.read_file(admin0_path)
    admin1 = gpd.read_file(admin1_path)
    
    # Filter China
    china = admin0[admin0["ADMIN"] == "China"]
    if "admin" in admin1.columns:
        provinces = admin1[admin1["admin"] == "China"]
    elif "ADM0NAME" in admin1.columns:
        provinces = admin1[admin1["ADM0NAME"] == "China"]
    else:
        provinces = gpd.GeoDataFrame(columns=admin1.columns, geometry=admin1.geometry.name)
    
    # Reproject to Asia Lambert Conformal Conic
    target_crs = "ESRI:102012"
    china = china.to_crs(target_crs)
    provinces = provinces.to_crs(target_crs)
    dsa_gdf = dsa_gdf.to_crs(target_crs)
    eor_gdf = eor_gdf.to_crs(target_crs)
    gas_gdf = gas_gdf.to_crs(target_crs)
    
    # Compute 250km buffers around gas fields
    buffer_250km = gas_gdf.buffer(250_000)
    buffer_union = unary_union(buffer_250km)
    
    # Find overlap with storage potential
    overlap_polygons = []
    for pt in pd.concat([dsa_gdf, eor_gdf]).geometry:
        if buffer_union.contains(pt):
            overlap_polygons.append(buffer_union)
    
    overlap_gdf = gpd.GeoDataFrame(
        geometry=[unary_union(overlap_polygons)],
        crs=target_crs
    ) if overlap_polygons else gpd.GeoDataFrame(geometry=[], crs=target_crs)
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Plot base layers
    china.boundary.plot(ax=ax, color="black", linewidth=1)
    provinces.boundary.plot(ax=ax, color="gray", linewidth=0.5, alpha=0.7)
    
    # Plot storage potential
    dsa_gdf.plot(ax=ax, color="green", markersize=25, alpha=0.7, label="DSA storage potential")
    eor_gdf.plot(ax=ax, color="orange", markersize=25, alpha=0.7, label="EOR storage potential")
    
    # Plot gas fields (size scaled by area if available)
    if "Shape_Area" in gas_gdf.columns:
        sizes = (gas_gdf["Shape_Area"] / gas_gdf["Shape_Area"].max()) * 200
    else:
        sizes = 50
    gas_gdf.plot(ax=ax, color="red", markersize=sizes, alpha=0.6, label="Gas fields")
    
    # Plot overlap region
    if not overlap_gdf.empty:
        overlap_gdf.plot(
            ax=ax, 
            facecolor="none", 
            edgecolor="blue",
            hatch="///", 
            linewidth=0, 
            alpha=0.5
        )
    
    # Create legend
    legend_elements = [
        Patch(facecolor="green", edgecolor="black", label="DSA storage potential"),
        Patch(facecolor="orange", edgecolor="black", label="EOR storage potential"),
        Patch(facecolor="red", edgecolor="black", label="Gas fields"),
        Patch(facecolor="none", edgecolor="blue", hatch="///", 
              label="Distance between gas fields and storage potential < 250km"),
    ]
    labels = ["\n".join(wrap(elem.get_label(), 30)) for elem in legend_elements]
    ax.legend(
        handles=legend_elements,
        labels=labels,
        loc="lower left",
        fontsize=9,
        handlelength=2,
        handleheight=1.5,
        labelspacing=1.2,
        frameon=True
    )

    ax.axis("off")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Map saved to {output_path}")


# Run function
plot_storage_potential_and_gas_fields(
    dsa_pipeline="./data/DSA-Pipeline.xlsx",
    eor_pipeline="./data/EOR_Pipeline.xlsx",
    china_oil_gas_csv="./data/china_gas_fields.csv",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="storage_gas_map.png"
)
