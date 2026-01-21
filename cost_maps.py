"""
Cost Maps Module
================
Functions for plotting China cost maps from Excel data files.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def plot_cost_maps(
    excel1_path, 
    excel2_path, 
    ne_admin0_dir, 
    ne_admin1_dir, 
    output_path="china_cost_maps.png",
    cost_type="avg"
):
    """
    Plot cost maps for 2025-2050 from two Excel files.

    Parameters:
        excel1_path: Path to first Excel file (e.g., DSA data)
        excel2_path: Path to second Excel file (e.g., EOR data)
        ne_admin0_dir: Path to Natural Earth admin0 shapefile directory
        ne_admin1_dir: Path to Natural Earth admin1 shapefile directory
        output_path: Output image path
        cost_type: "avg", "min", or "max" to select cost column type

    Color coding:
        - RED: Negative costs (clipped to 0)
        - GREY: Zero cost with zero injection rate but positive storage potential
        - Viridis colormap: Positive costs (converted CNY to USD)
    """

    # Percentile threshold for color scale normalization
    year_percentiles = {
        2025: 0.90,
        2030: 0.90,
        2035: 0.90,
        2040: 0.90,
        2045: 0.90,
        2050: 0.90
    }

    # Load Excel files
    df1 = pd.read_excel(excel1_path)
    df2 = pd.read_excel(excel2_path)

    # Validate required columns
    for dfi, tag in [(df1, "file 1"), (df2, "file 2")]:
        for req in ["X", "Y", "Injection Rate", "Storage Potential"]:
            if req not in dfi.columns:
                raise KeyError(f"Missing column '{req}' in {tag}")

    # Check for overlapping points
    overlap = pd.merge(df1[['X','Y']], df2[['X','Y']], on=['X','Y'])
    if not overlap.empty:
        print(f"Warning: {len(overlap)} clashing points found (present in both files)")
    else:
        print("No clashing points found between the two files.")

    # Merge datasets on coordinates
    cost_years = [2025, 2030, 2035, 2040, 2045, 2050]
    merged = pd.merge(df1, df2, on=["X", "Y"], how="outer", suffixes=("_1", "_2"))

    # Merge Injection Rate and Storage Potential for grey point classification
    merged["Injection Rate_merged"] = merged[["Injection Rate_1", "Injection Rate_2"]].min(axis=1, skipna=True)
    merged["Storage Potential_merged"] = merged[["Storage Potential_1", "Storage Potential_2"]].max(axis=1, skipna=True)

    # Compute per-year costs (use minimum when both sources have data)
    for year in cost_years:
        col1 = f"Cost_{cost_type}_{year}_1"
        col2 = f"Cost_{cost_type}_{year}_2"
        if col1 not in merged.columns and col2 not in merged.columns:
            raise KeyError(f"Neither '{col1}' nor '{col2}' present in inputs.")
        merged[f"Cost_{cost_type}_{year}_orig"] = merged[[c for c in [col1, col2] if c in merged.columns]].min(axis=1, skipna=True)
        merged[f"Cost_{cost_type}_{year}"] = merged[f"Cost_{cost_type}_{year}_orig"].clip(lower=0)

    # Convert CNY to USD
    CNY_to_USD = 1 / 7.183
    for year in cost_years:
        merged[f"Cost_{cost_type}_{year}_USD"] = merged[f"Cost_{cost_type}_{year}"] * CNY_to_USD

    # Load shapefiles
    admin0 = gpd.read_file(f"{ne_admin0_dir}/ne_110m_admin_0_countries.shp")
    china_outline = admin0[admin0["NAME"] == "China"].copy()
    if china_outline.empty:
        raise ValueError("China outline not found in admin0 shapefile.")

    admin1 = gpd.read_file(f"{ne_admin1_dir}/ne_50m_admin_1_states_provinces.shp")
    if "admin" in admin1.columns:
        china_provinces = admin1[admin1["admin"] == "China"].copy()
    elif "ADM0NAME" in admin1.columns:
        china_provinces = admin1[admin1["ADM0NAME"] == "China"].copy()
    else:
        china_provinces = gpd.GeoDataFrame(columns=admin1.columns, geometry=admin1.geometry.name)

    # Clean and reproject geometries
    china_provinces = china_provinces[china_provinces.is_valid & ~china_provinces.geometry.is_empty]
    if not china_provinces.empty:
        china_provinces = china_provinces.explode(index_parts=False)
        china_provinces["geometry"] = china_provinces.geometry.intersection(china_outline.unary_union)

    china_outline = china_outline.to_crs(epsg=3857)
    if not china_provinces.empty:
        china_provinces = china_provinces.to_crs(epsg=3857)

    # Create GeoDataFrame of points
    gdf_points = gpd.GeoDataFrame(
        merged,
        geometry=gpd.points_from_xy(merged["X"], merged["Y"]),
        crs="EPSG:4326"
    ).to_crs(epsg=3857)

    # Create figure with 2x3 subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), constrained_layout=True)
    axes = axes.ravel()

    def get_subplot_label(year):
        """Return subplot label (a)-(f) based on year."""
        labels = {2025: "(a)", 2030: "(b)", 2035: "(c)", 2040: "(d)", 2045: "(e)", 2050: "(f)"}
        return labels.get(year, "")

    for i, year in enumerate(cost_years):
        col_usd = f"Cost_{cost_type}_{year}_USD"
        col_orig = f"Cost_{cost_type}_{year}_orig"

        # Define point masks
        mask_red = gdf_points[col_orig] < 0
        mask_grey = (
            (gdf_points[col_orig] == 0) & 
            (gdf_points["Injection Rate_merged"].fillna(float("inf")) == 0) & 
            (gdf_points["Storage Potential_merged"].fillna(0) > 0)
        )
        mask_pos = gdf_points[col_usd] > 0

        # Plot boundaries
        china_outline.boundary.plot(ax=axes[i], linewidth=1.5, color="black")
        if not china_provinces.empty:
            china_provinces.plot(ax=axes[i], facecolor="none", edgecolor="gray", linewidth=1)

        # Plot positive cost points with colormap
        if not gdf_points[mask_pos].empty:
            percentile = year_percentiles[year]
            vals = gdf_points.loc[mask_pos, col_usd]
            vmin = 0
            vmax = vals.quantile(percentile)
            
            plot_vals = vals.clip(upper=vmax)
            gdf_plot = gdf_points[mask_pos].copy()
            gdf_plot[f"{col_usd}_plot"] = plot_vals
            
            gdf_plot.plot(
                ax=axes[i],
                column=f"{col_usd}_plot",
                cmap="viridis_r",
                markersize=15,
                legend=False,
                vmin=vmin,
                vmax=vmax
            )
            
            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=plt.Normalize(vmin=vmin, vmax=vmax))
            sm._A = []
            cbar = fig.colorbar(sm, ax=axes[i], shrink=0.8)
            cbar.locator = MaxNLocator(nbins=6, prune=None)
            cbar.update_ticks()
            ticks = cbar.get_ticks()
            tick_labels = [f"{t:.0f}" for t in ticks[:-1]] + [f">={int(vmax)}"]
            cbar.set_ticklabels(tick_labels)
            cbar.set_label("Cost (USD/km)")

        # Plot negative cost points (red)
        if not gdf_points[mask_red].empty:
            gdf_points[mask_red].plot(ax=axes[i], color="red", markersize=15)

        # Plot special zero points (grey)
        if not gdf_points[mask_grey].empty:
            gdf_points[mask_grey].plot(ax=axes[i], color="grey", markersize=15)

        # Format subplot
        bounds = china_outline.total_bounds
        axes[i].set_xlim(bounds[0], bounds[2])
        axes[i].set_ylim(bounds[1], bounds[3])
        axes[i].set_title(get_subplot_label(year))
        axes[i].set_aspect('equal')
        axes[i].set_xticks([])
        axes[i].set_yticks([])

    plt.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"[OK] Plot saved to {output_path}")
