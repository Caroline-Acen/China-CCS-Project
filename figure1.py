"""
Figure 1: Injection Maps
Plots China maps showing Storage Potential and Injection Rate from 4 Excel files.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def plot_injection_maps(
    excel1, excel2, excel3, excel4,
    admin0_path, admin1_path,
    output_path="china_injection_maps.png"
):
    """
    Plot China maps from 4 Excel files showing Storage Potential and Injection Rate.

    Parameters:
        excel1-4: Paths to Excel files (DSA/EOR Pipeline/Tank data)
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path

    Color coding:
        - RED: Negative values (clipped to 0) or true zeros
        - GREY: Zero value with zero injection but positive storage potential
        - Viridis colormap: Positive values
    """

    # Percentile threshold for color scale
    subplot_percentiles = {0: 0.95, 1: 0.95, 2: 0.95, 3: 0.95}

    # Load country outline
    admin0 = gpd.read_file(admin0_path)
    china_outline = admin0[admin0["NAME"] == "China"].copy()
    if china_outline.empty:
        raise ValueError("China outline not found in admin0 shapefile.")

    # Load provincial boundaries
    admin1 = gpd.read_file(admin1_path)
    if "admin" in admin1.columns:
        china_provinces = admin1[admin1["admin"] == "China"].copy()
    elif "ADM0NAME" in admin1.columns:
        china_provinces = admin1[admin1["ADM0NAME"] == "China"].copy()
    else:
        print("Warning: Could not find suitable column for China provinces.")
        china_provinces = gpd.GeoDataFrame(columns=admin1.columns, geometry=admin1.geometry.name)

    # Clean geometries
    china_provinces = china_provinces[china_provinces.is_valid & ~china_provinces.geometry.is_empty]
    if not china_provinces.empty:
        china_provinces = china_provinces.explode(index_parts=False)

    # Reproject to Web Mercator
    china_outline = china_outline.to_crs(epsg=3857)
    if not china_provinces.empty:
        china_provinces = china_provinces.to_crs(epsg=3857)

    # Read Excel files into GeoDataFrames
    gdfs = []
    for path in [excel1, excel2, excel3, excel4]:
        df = pd.read_excel(path)

        # Keep original values for mask classification
        df["Storage Potential_orig"] = df["Storage Potential"]
        df["Injection Rate_orig"] = df["Injection Rate"]

        # Clip negative values to zero
        df["Storage Potential"] = df["Storage Potential"].clip(lower=0)
        df["Injection Rate"] = df["Injection Rate"].clip(lower=0)

        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(df["X"], df["Y"]),
            crs="EPSG:4326"
        )
        gdfs.append(gdf.to_crs(epsg=3857))

    # Columns to plot for each subplot
    columns = ["Storage Potential", "Storage Potential", "Injection Rate", "Injection Rate"]

    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(18, 14), constrained_layout=True)
    axes = axes.ravel()
    titles = ["(a)", "(b)", "(c)", "(d)"]

    for i, (ax, gdf, title, col) in enumerate(zip(axes, gdfs, titles, columns)):
        label = "Storage Potential (Mt)" if col == "Storage Potential" else "Injection Rate (Mt/year)"
        orig_col = f"{col}_orig"

        # Define masks
        mask_pos = gdf[col] > 0
        mask_red = gdf[orig_col] <= 0
        mask_grey = (
            (gdf[orig_col] == 0) &
            (gdf["Injection Rate_orig"] == 0) &
            (gdf["Storage Potential_orig"] != 0)
        )
        mask_red = mask_red & ~mask_grey

        # Plot boundaries
        china_outline.boundary.plot(ax=ax, linewidth=1.5, color="black")
        if not china_provinces.empty:
            china_provinces.plot(ax=ax, facecolor="none", edgecolor="gray", linewidth=1)

        # Plot positive value points
        if not gdf[mask_pos].empty:
            percentile = subplot_percentiles[i]
            vals = gdf[mask_pos][col]
            vmin = 0
            vmax = vals.quantile(percentile)
            
            plot_vals = vals.clip(upper=vmax)
            gdf_plot = gdf[mask_pos].copy()
            gdf_plot[f"{col}_plot"] = plot_vals
            
            gdf_plot.plot(
                ax=ax,
                column=f"{col}_plot",
                cmap="viridis_r",
                markersize=15,
                legend=False,
                vmin=vmin,
                vmax=vmax
            )

            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=plt.Normalize(vmin=vmin, vmax=vmax))
            sm._A = []
            cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.7)
            cbar.locator = MaxNLocator(nbins=6, prune=None)
            cbar.update_ticks()
            ticks = cbar.get_ticks()
            tick_labels = [f"{t:.0f}" for t in ticks[:-1]] + [f">={int(vmax)}"]
            cbar.set_ticklabels(tick_labels)
            cbar.set_label(label)

        # Plot negative/zero points (red)
        if not gdf[mask_red].empty:
            gdf[mask_red].plot(ax=ax, color="red", markersize=15)

        # Plot special zero points (grey)
        if not gdf[mask_grey].empty:
            gdf[mask_grey].plot(ax=ax, color="grey", markersize=15)

        # Format subplot
        bounds = china_outline.total_bounds
        ax.set_xlim(bounds[0], bounds[2])
        ax.set_ylim(bounds[1], bounds[3])
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    plt.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"[OK] Plot saved to {output_path}")


# Run function
plot_injection_maps(
    excel1="./data/DSA-Pipeline.xlsx",
    excel2="./data/EOR_Pipeline.xlsx",
    excel3="./data/DSA-Tank.xlsx",
    excel4="./data/EOR_Tank.xlsx",
    admin0_path="./data/natural_earth/admin0",
    admin1_path="./data/natural_earth/admin1",
    output_path="injection_maps.png"
)
