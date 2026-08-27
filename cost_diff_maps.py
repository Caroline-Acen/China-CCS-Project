"""
Cost Difference Maps Module
============================
Functions for plotting China cost difference maps comparing Pipeline vs Tank transport.
Positive values (blue) = Pipeline cheaper, Negative values (red) = Tank cheaper
Uses percentage difference: ((Tank - Pipeline) / Tank) * 100
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from config import CNY_TO_USD, manuscript_font_bundle


def _classify_cost_points(cost_series, injection_series, storage_series):
    """Classify points into positive, exhausted, and undetermined buckets for masks.

    Matches cost_maps._classify_cost_points: red when cost ≤ 0 / missing / NaN /
    #DIV/0!; deep green when injection = 0.
    """
    _ = storage_series
    injection_series = injection_series.fillna(0)

    cost_numeric = pd.to_numeric(cost_series, errors="coerce")
    mask_cost_bad = cost_numeric.le(0) | cost_numeric.isna()
    mask_grey = injection_series.eq(0)
    mask_red = mask_cost_bad & ~mask_grey
    mask_pos = cost_numeric.gt(0) & ~mask_grey
    return mask_pos, mask_red, mask_grey


DIFF_CMAP = LinearSegmentedColormap.from_list(
    "tank_warm_pipeline_cool",
    ["#5a1d11", "#d8e68b", "#f7f7f7", "#90d6d0", "#2b2c8e"],
)


def plot_cost_difference_maps(
    pipeline_excel1_path,
    pipeline_excel2_path,
    tank_excel1_path, 
    tank_excel2_path,
    ne_admin0_dir, 
    ne_admin1_dir, 
    output_path="china_cost_difference_maps.png",
    cost_type="avg",
    colorbar=1,
    apply_cny_conversion=False,
):
    """
    Plot cost difference maps ((Tank - Pipeline) / Tank) * 100 for 2025-2050 from Excel files.
    
    Positive difference (blue/cyan): Pipeline transport is more economical
    Negative difference (yellow/brown): Tank transport is more economical

    Parameters:
        pipeline_excel1_path: Path to first pipeline Excel file (e.g., DSA-Pipeline)
        pipeline_excel2_path: Path to second pipeline Excel file (e.g., EOR-Pipeline)
        tank_excel1_path: Path to first tank Excel file (e.g., DSA-Tank)
        tank_excel2_path: Path to second tank Excel file (e.g., EOR-Tank)
        ne_admin0_dir: Path to Natural Earth admin0 shapefile directory
        ne_admin1_dir: Path to Natural Earth admin1 shapefile directory
        output_path: Output image path
        cost_type: "avg", "min", or "max" to select cost column type
        colorbar: 1 (default) for one unified colorbar, 2 for two row-wise colorbars

    Color coding:
        - RED: Exhausted storage potential (cost ≤ 0, missing/NaN, or #DIV/0!)
        - BLUE shades: Pipeline transport is more economical (positive % difference)
        - WHITE: Similar costs (near zero % difference)
        - BROWN shades: Tank transport is more economical (negative % difference)
    """

    # Percentile threshold for color scale normalization
    percentile = 0.95

    # Load Pipeline Excel files
    df_pipe1 = pd.read_excel(pipeline_excel1_path)
    df_pipe2 = pd.read_excel(pipeline_excel2_path)

    # Load Tank Excel files
    df_tank1 = pd.read_excel(tank_excel1_path)
    df_tank2 = pd.read_excel(tank_excel2_path)

    # Validate required columns
    for dfi, tag in [(df_pipe1, "pipeline file 1"), (df_pipe2, "pipeline file 2"),
                     (df_tank1, "tank file 1"), (df_tank2, "tank file 2")]:
        for req in ["X", "Y", "Injection Rate", "Storage Potential"]:
            if req not in dfi.columns:
                raise KeyError(f"Missing column '{req}' in {tag}")

    # Merge pipeline datasets on coordinates
    cost_years = [2025, 2030, 2035, 2040, 2045, 2050]
    merged_pipe = pd.merge(df_pipe1, df_pipe2, on=["X", "Y"], how="outer", suffixes=("_1", "_2"))

    # Merge tank datasets on coordinates
    merged_tank = pd.merge(df_tank1, df_tank2, on=["X", "Y"], how="outer", suffixes=("_1", "_2"))

    # Merge Injection Rate and Storage Potential for exhausted-point classification
    merged_pipe["Injection Rate_pipe"] = merged_pipe[["Injection Rate_1", "Injection Rate_2"]].min(axis=1, skipna=True)
    merged_pipe["Storage Potential_pipe"] = merged_pipe[["Storage Potential_1", "Storage Potential_2"]].max(axis=1, skipna=True)
    merged_tank["Injection Rate_tank"] = merged_tank[["Injection Rate_1", "Injection Rate_2"]].min(axis=1, skipna=True)
    merged_tank["Storage Potential_tank"] = merged_tank[["Storage Potential_1", "Storage Potential_2"]].max(axis=1, skipna=True)

    # Compute per-year costs for pipeline (use minimum when both sources have data)
    for year in cost_years:
        col1 = f"Cost_{cost_type}_{year}_1"
        col2 = f"Cost_{cost_type}_{year}_2"
        if col1 in merged_pipe.columns or col2 in merged_pipe.columns:
            merged_pipe[f"Cost_{cost_type}_{year}_pipe_orig"] = merged_pipe[[c for c in [col1, col2] if c in merged_pipe.columns]].min(axis=1, skipna=True)

    # Compute per-year costs for tank
    for year in cost_years:
        col1 = f"Cost_{cost_type}_{year}_1"
        col2 = f"Cost_{cost_type}_{year}_2"
        if col1 in merged_tank.columns or col2 in merged_tank.columns:
            merged_tank[f"Cost_{cost_type}_{year}_tank_orig"] = merged_tank[[c for c in [col1, col2] if c in merged_tank.columns]].min(axis=1, skipna=True)

    # Merge pipeline and tank data on coordinates
    merged = pd.merge(merged_pipe[
                          ["X", "Y", "Injection Rate_pipe", "Storage Potential_pipe"]
                          + [f"Cost_{cost_type}_{year}_pipe_orig" for year in cost_years]
                      ], 
                      merged_tank[
                          ["X", "Y", "Injection Rate_tank", "Storage Potential_tank"]
                          + [f"Cost_{cost_type}_{year}_tank_orig" for year in cost_years]
                      ], 
                      on=["X", "Y"], how="outer")

    # Optional CNY → USD conversion (legacy data/ files only)
    fx = CNY_TO_USD if apply_cny_conversion else 1.0

    # Calculate cost differences on valid points only:
    # - exclude exhausted points (non-positive cost, except undetermined zero-injection storage cells)
    # - include undetermined points in the valid pool when both costs are positive
    # - exclude non-positive costs from difference calculation
    # Formula: Cost_diff = ((Cost_tank_USD - Cost_pipe_USD) / Cost_tank_USD) * 100
    for year in cost_years:
        pipe_col_orig = f"Cost_{cost_type}_{year}_pipe_orig"
        tank_col_orig = f"Cost_{cost_type}_{year}_tank_orig"

        merged[f"Cost_{cost_type}_{year}_pipe_USD"] = merged[pipe_col_orig].clip(lower=0) * fx
        merged[f"Cost_{cost_type}_{year}_tank_USD"] = merged[tank_col_orig].clip(lower=0) * fx

        pipe_pos, pipe_red, pipe_grey = _classify_cost_points(
            merged[pipe_col_orig],
            merged["Injection Rate_pipe"],
            merged["Storage Potential_pipe"],
        )
        tank_pos, tank_red, tank_grey = _classify_cost_points(
            merged[tank_col_orig],
            merged["Injection Rate_tank"],
            merged["Storage Potential_tank"],
        )
        mask_grey_special = pipe_grey | tank_grey
        mask_red = (pipe_red | tank_red) & (~mask_grey_special)

        mask_valid = (~mask_red) & pipe_pos & tank_pos
        mask_outline = (~mask_valid) & (~mask_red)

        merged[f"Mask_red_{year}"] = mask_red
        merged[f"Mask_outline_{year}"] = mask_outline

        merged[f"Cost_diff_{year}_pct"] = np.where(
            mask_valid,
            ((merged[f"Cost_{cost_type}_{year}_tank_USD"] - merged[f"Cost_{cost_type}_{year}_pipe_USD"]) / merged[f"Cost_{cost_type}_{year}_tank_USD"]) * 100,
            np.nan
        )

    # Calculate global scale across all years for diverging colormap
    all_diffs = []
    for year in cost_years:
        col_diff = f"Cost_diff_{year}_pct"
        valid_diffs = merged[col_diff].dropna()
        # Filter out extreme outliers for better visualization
        valid_diffs = valid_diffs[np.isfinite(valid_diffs)]
        all_diffs.extend(valid_diffs.tolist())
    
    if all_diffs:
        # Use symmetric scale around zero
        abs_max = np.percentile(np.abs(all_diffs), percentile * 100)
        global_vmin = -abs_max
        global_vmax = abs_max
    else:
        global_vmin = -100
        global_vmax = 100

    # If user requests two colorbars, compute a symmetric scale for each row separately
    if colorbar == 2:
        row1_vals = []
        row2_vals = []
        for year in cost_years[:3]:
            col_diff = f"Cost_diff_{year}_pct"
            vals = merged[col_diff].dropna()
            vals = vals[np.isfinite(vals)]
            row1_vals.extend(vals.tolist())
        for year in cost_years[3:]:
            col_diff = f"Cost_diff_{year}_pct"
            vals = merged[col_diff].dropna()
            vals = vals[np.isfinite(vals)]
            row2_vals.extend(vals.tolist())

        if row1_vals:
            row1_abs_max = np.percentile(np.abs(row1_vals), percentile * 100)
            row1_vmin, row1_vmax = -row1_abs_max, row1_abs_max
        else:
            row1_vmin, row1_vmax = -100, 100

        if row2_vals:
            row2_abs_max = np.percentile(np.abs(row2_vals), percentile * 100)
            row2_vmin, row2_vmax = -row2_abs_max, row2_abs_max
        else:
            row2_vmin, row2_vmax = -100, 100
    else:
        row1_vmin = row2_vmin = global_vmin
        row1_vmax = row2_vmax = global_vmax

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

    # Precompute outline-only locations across years and restrict to the
    # two geographic regions that should retain black boundaries.
    mask_cols = [f"Mask_outline_{y}" for y in cost_years if f"Mask_outline_{y}" in merged.columns]
    if mask_cols:
        merged["Mask_outline_any"] = merged[mask_cols].any(axis=1)
    else:
        merged["Mask_outline_any"] = False

    # Only keep outline points in the two requested map regions.
    outline_box1 = (
        (merged["X"] >= 108) &
        (merged["X"] <= 118) &
        (merged["Y"] >= 15) &
        (merged["Y"] <= 24)
    )
    outline_box2 = (
        (merged["X"] >= 120) &
        (merged["X"] <= 129) &
        (merged["Y"] >= 25) &
        (merged["Y"] <= 34)
    )
    merged["Mask_outline_any"] = merged["Mask_outline_any"] & (outline_box1 | outline_box2)

    # gdf for the 'any-year' outline-only locations
    gdf_outline_any = gdf_points[merged["Mask_outline_any"].values]
    # Build a single buffered outline from excluded points so we can draw
    # a consistent hollow shape on every subplot. Use a modest buffer so the
    # boundary follows cluster shapes more closely.
    outline_buffer_m = 25000  # 25 km buffer to create contiguous outlines
    if not gdf_outline_any.empty:
        try:
            union_geom = gdf_outline_any.geometry.unary_union
            outline_geom = union_geom.buffer(outline_buffer_m).boundary
            gdf_outline = gpd.GeoDataFrame(geometry=[outline_geom], crs=gdf_outline_any.crs)
        except Exception:
            gdf_outline = gdf_outline_any
    else:
        gdf_outline = gpd.GeoDataFrame(geometry=[], crs=gdf_points.crs)

    # Create figure with 2x3 subplots and space for legend
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    mf = manuscript_font_bundle(fig.get_figwidth())
    axes = axes.ravel()

    def get_subplot_label(year):
        """Return subplot label (a)-(f) based on year."""
        labels = {2025: "(a)", 2030: "(b)", 2035: "(c)", 2040: "(d)", 2045: "(e)", 2050: "(f)"}
        return labels.get(year, "")

    for i, year in enumerate(cost_years):
        col_diff = f"Cost_diff_{year}_pct"

        # Define point mask for valid data
        mask_valid = ~gdf_points[col_diff].isna() & np.isfinite(gdf_points[col_diff])
        mask_red = gdf_points[f"Mask_red_{year}"] == True
        # Exclude any-location that was outline-only in any year from being filled
        if "Mask_outline_any" in merged.columns:
            not_outline_any = ~merged["Mask_outline_any"].values
            mask_valid = mask_valid & not_outline_any
            mask_red = mask_red & not_outline_any

        # Plot boundaries
        china_outline.boundary.plot(ax=axes[i], linewidth=1.5, color="black")
        if not china_provinces.empty:
            china_provinces.plot(ax=axes[i], facecolor="none", edgecolor="gray", linewidth=1)

        # Plot cost difference points with diverging colormap
        if not gdf_points[mask_valid].empty:
            if colorbar == 2:
                vmin_use = row1_vmin if i < 3 else row2_vmin
                vmax_use = row1_vmax if i < 3 else row2_vmax
            else:
                vmin_use = global_vmin
                vmax_use = global_vmax

            gdf_points[mask_valid].plot(
                ax=axes[i],
                column=col_diff,
                cmap=DIFF_CMAP,
                markersize=15,
                legend=False,
                vmin=vmin_use,
                vmax=vmax_use
            )

        # Plot exhausted/non-positive points (red category)
        if not gdf_points[mask_red].empty:
            gdf_points[mask_red].plot(ax=axes[i], color="red", markersize=15)

        # Plot undetermined (mask_grey) as hollow markers with black edges (no fill)
        # Draw the undetermined outline polygon (no fill) on each subplot so
        # the region is visible but not filled. If we have an outline polygon,
        # plot it; otherwise fall back to point markers.
        if not gdf_outline.empty:
            gdf_outline.plot(ax=axes[i], facecolor="none", edgecolor="black", linewidth=1.5, zorder=13)
        elif not gdf_outline_any.empty:
            gdf_outline_any.plot(
                ax=axes[i],
                marker="o",
                facecolor="none",
                edgecolor="black",
                linewidth=1.0,
                markersize=50,
                zorder=13,
            )

        # Re-plot boundaries on top so empty areas keep their outlines even when no points are drawn
        china_outline.boundary.plot(ax=axes[i], linewidth=1.5, color="black", zorder=12)
        if not china_provinces.empty:
            china_provinces.plot(ax=axes[i], facecolor="none", edgecolor="gray", linewidth=1, zorder=11)

        # Format subplot
        bounds = china_outline.total_bounds
        axes[i].set_xlim(bounds[0], bounds[2])
        axes[i].set_ylim(bounds[1], bounds[3])
        axes[i].set_title(get_subplot_label(year), fontsize=mf["title"], fontweight="bold")
        axes[i].set_aspect('equal')
        axes[i].set_xticks([])
        axes[i].set_yticks([])

    # Add colorbar(s) and legend
    # Tighten left margin so the map grid doesn't leave a large blank strip.
    plt.subplots_adjust(left=0.04, right=0.88, hspace=0.15, wspace=0.1)

    if colorbar == 1:
        # Unified colorbar for all panels
        cbar_ax = fig.add_axes([0.91, 0.15, 0.015, 0.7])  # [left, bottom, width, height]
        sm = plt.cm.ScalarMappable(cmap=DIFF_CMAP, norm=plt.Normalize(vmin=global_vmin, vmax=global_vmax))
        sm._A = []
        cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical")
        cbar.locator = MaxNLocator(nbins=8)
        cbar.update_ticks()
        cbar.set_label("Cost Difference: ((Tank − Pipeline) / Tank) × 100 (%)", fontsize=mf["colorbar_label"])
        cbar.ax.tick_params(labelsize=mf["colorbar_tick"])
    else:
        # Two colorbars: one for top row, one for bottom row
        sm1 = plt.cm.ScalarMappable(cmap=DIFF_CMAP, norm=plt.Normalize(vmin=row1_vmin, vmax=row1_vmax))
        sm1._A = []
        cbar_ax1 = fig.add_axes([0.91, 0.53, 0.015, 0.35])
        cbar1 = fig.colorbar(sm1, cax=cbar_ax1, orientation="vertical")
        cbar1.locator = MaxNLocator(nbins=8)
        cbar1.update_ticks()
        cbar1.ax.tick_params(labelsize=mf["colorbar_tick"])

        sm2 = plt.cm.ScalarMappable(cmap=DIFF_CMAP, norm=plt.Normalize(vmin=row2_vmin, vmax=row2_vmax))
        sm2._A = []
        cbar_ax2 = fig.add_axes([0.91, 0.12, 0.015, 0.35])
        cbar2 = fig.colorbar(sm2, cax=cbar_ax2, orientation="vertical")
        cbar2.locator = MaxNLocator(nbins=8)
        cbar2.update_ticks()
        cbar2.ax.tick_params(labelsize=mf["colorbar_tick"])

        fig.text(
            0.975,
            0.50,
            "Cost Difference: ((Tank − Pipeline) / Tank) × 100 (%)",
            rotation=90,
            va="center",
            ha="center",
            fontsize=mf["colorbar_label"],
        )
    
    # Add legend with adaptive wrapping to avoid clipping at figure edges
    legend_labels = [
         'RED: Exhausted storage potential',
           'BLUE: Pipeline transport is the more cost-effective option',
            'BROWN: Tank transport is the more cost-effective option'
    ]

    legend_handles = [
         Line2D([0], [0], marker='o', color='none', markerfacecolor='red', markeredgecolor='red', markersize=mf["legend_marker"],
               label=legend_labels[0]),
        Line2D([0], [0], marker='o', color='none', markerfacecolor='blue', markeredgecolor='blue', markersize=mf["legend_marker"],
             label=legend_labels[1]),
         Line2D([0], [0], marker='o', color='none', markerfacecolor='#a0782b', markeredgecolor='#a0782b', markersize=mf["legend_marker"],
             label=legend_labels[2])
    ]

    # Estimate how many legend entries fit on one row and wrap when needed
    figure_width_px = fig.get_size_inches()[0] * fig.dpi
    usable_width_px = figure_width_px * 0.86
    item_widths = [70 + (len(label) * 8.5) for label in legend_labels]

    legend_ncol = 1
    for candidate_ncol in range(len(legend_labels), 0, -1):
        row_widths = []
        for row_start in range(0, len(legend_labels), candidate_ncol):
            row_widths.append(sum(item_widths[row_start:row_start + candidate_ncol]))
        if max(row_widths) <= usable_width_px:
            legend_ncol = candidate_ncol
            break

    legend_rows = int(np.ceil(len(legend_labels) / legend_ncol))
    legend_height = 0.08 if legend_rows == 1 else (0.05 + 0.035 * legend_rows)

    legend_ax = fig.add_axes([0.02, 0.01, 0.86, legend_height])
    legend_ax.axis('off')

    legend_ax.legend(
        handles=legend_handles,
        loc='center',
        ncol=legend_ncol,
        frameon=False,
        fontsize=mf["legend"],
        handletextpad=0.5,
        columnspacing=1.4
    )

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[OK] Cost difference map saved to {output_path}")
