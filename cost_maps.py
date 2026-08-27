"""
Cost Maps Module
================
Functions for plotting China cost maps from Excel data files.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

from config import (
    CNY_TO_USD,
    COST_SCENARIOS,
    COST_YEARS,
    manuscript_font_bundle,
)


# Reference-style sequential palette: deep blue -> cyan -> yellow -> brown
# for progressively higher positive costs.
COST_CMAP = LinearSegmentedColormap.from_list(
    "cost_blue_cyan_yellow_brown",
    ["#2b2c8e", "#5b88c5", "#90d6d0", "#d8e68b", "#a0782b", "#5a1d11"],
)

# Distinct category colors
COLOR_EXHAUSTED = "red"
COLOR_UNDETERMINED = "#00441b"  # deep green


def _classify_cost_points(cost_series, injection_series, storage_series):
    """Classify points into positive, exhausted, and undetermined buckets.

    Red (exhausted storage potential) when cost is ≤ 0, missing/NaN, or
    non-numeric (e.g. Excel #DIV/0!), except undetermined zero-injection cells.

    Deep green (undetermined): injection = 0
      (e.g. baseline 2025 zero injection).

    Colormap: cost > 0 (and not red / undetermined).

    storage_series is kept for call-site compatibility; exhaustion is inferred
    from the cost column (static Storage Potential is rarely exactly 0).
    """
    _ = storage_series
    injection_series = injection_series.fillna(0)

    # #DIV/0! and other Excel errors become NaN via coerce.
    cost_numeric = pd.to_numeric(cost_series, errors="coerce")
    mask_cost_bad = cost_numeric.le(0) | cost_numeric.isna()

    mask_grey = injection_series.eq(0)
    mask_red = mask_cost_bad & ~mask_grey
    mask_pos = cost_numeric.gt(0) & ~mask_grey
    return mask_pos, mask_red, mask_grey


def compute_cost_statistics(
    excel_path,
    apply_cny_conversion=False,
    output_csv=None,
    label=None,
):
    """
    Compute minimum, mean, maximum, and percentile statistics for each
    year and cost scenario (min / avg / max).

    Returns a DataFrame formatted like Supplementary Table 1.
    """
    df = pd.read_excel(excel_path)
    rows = []

    scenario_titles = {"min": "Minimum", "avg": "Average", "max": "Maximum"}
    stat_funcs = {
        "Minimum": lambda s: s.min(),
        "Mean": lambda s: s.mean(),
        "Maximum": lambda s: s.max(),
        "5th percentile": lambda s: s.quantile(0.05),
        "10th percentile": lambda s: s.quantile(0.10),
        "95th percentile": lambda s: s.quantile(0.95),
    }

    for scenario in COST_SCENARIOS:
        for stat_name, stat_fn in stat_funcs.items():
            row = {
                "Scenario group": scenario_titles[scenario],
                "Statistic": stat_name,
            }
            for year in COST_YEARS:
                col = f"Cost_{scenario}_{year}"
                if col not in df.columns:
                    raise KeyError(f"Missing column '{col}' in {excel_path}")
                vals = df.loc[df[col] > 0, col].astype(float)
                if apply_cny_conversion:
                    vals = vals * CNY_TO_USD
                row[year] = stat_fn(vals) if not vals.empty else float("nan")
            rows.append(row)

    stats_df = pd.DataFrame(rows)
    if label:
        stats_df.insert(0, "Dataset", label)
    if output_csv:
        stats_df.to_csv(output_csv, index=False)
        print(f"[OK] Cost statistics saved to {output_csv}")
    return stats_df


def plot_cost_maps(
    excel1_path, 
    excel2_path, 
    ne_admin0_dir, 
    ne_admin1_dir, 
    output_path="china_cost_maps.png",
    cost_type="avg",
    colorbar=1,
    apply_cny_conversion=False,
    cost_unit_label="USD/tCO₂",
    stats_csv_path=None,
    stats_label=None,
    outline_coverage_paths=None,
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
        outline_coverage_paths: Optional Excel paths whose XY footprint defines
            persistent area outlines (even when plotted costs are empty/NaN).
            Defaults to excel1_path and excel2_path. For EOR-only maps, pass the
            companion DSA file so empty EOR regions keep the same boundaries.

    Color coding:
        - RED: Exhausted storage potential (cost ≤ 0, missing/NaN, or #DIV/0!)
        - DEEP GREEN: Undetermined injection rate (injection = 0)
        - Blue-cyan-yellow-brown colormap: Positive costs

    Parameter `colorbar`:
        - 1 (default): single unified colorbar for all panels (current behavior)
        - 2: one colorbar for the first row and another for the second row
    """

    # Percentile threshold for color scale normalization
    percentile = 0.90

    if stats_csv_path or stats_label:
        compute_cost_statistics(
            excel1_path,
            apply_cny_conversion=apply_cny_conversion,
            output_csv=stats_csv_path,
            label=stats_label or excel1_path,
        )

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

    # Merge Injection Rate and Storage Potential for undetermined point classification
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

    # Optional CNY → USD conversion (legacy data/ files only)
    for year in cost_years:
        usd_vals = merged[f"Cost_{cost_type}_{year}"]
        if apply_cny_conversion:
            usd_vals = usd_vals * CNY_TO_USD
        merged[f"Cost_{cost_type}_{year}_USD"] = usd_vals

    # Calculate global / per-row scale across years
    all_positive_costs = []
    for year in cost_years:
        col_usd = f"Cost_{cost_type}_{year}_USD"
        positive_vals = merged[merged[col_usd] > 0][col_usd]
        all_positive_costs.extend(positive_vals.tolist())

    if all_positive_costs:
        global_vmin = 0
        global_vmax = pd.Series(all_positive_costs).quantile(percentile)
    else:
        global_vmin = 0
        global_vmax = 1

    # If user requests two colorbars, compute a vmax for each row separately
    if colorbar == 2:
        # first row: years 0-2, second row: years 3-5
        row1_vals = []
        row2_vals = []
        for year in cost_years[:3]:
            col_usd = f"Cost_{cost_type}_{year}_USD"
            row1_vals.extend(merged[merged[col_usd] > 0][col_usd].tolist())
        for year in cost_years[3:]:
            col_usd = f"Cost_{cost_type}_{year}_USD"
            row2_vals.extend(merged[merged[col_usd] > 0][col_usd].tolist())

        if row1_vals:
            row1_vmax = pd.Series(row1_vals).quantile(percentile)
        else:
            row1_vmax = 1
        if row2_vals:
            row2_vmax = pd.Series(row2_vals).quantile(percentile)
        else:
            row2_vmax = 1
    else:
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

    # Build the outline-only region from undetermined (mask_grey) points,
    # missing-cost coverage, and any reference footprint points. Restrict to the
    # two target boxes so only those coastal/offshore regions retain the black
    # boundary treatment — including on EOR maps where those areas have no fill.
    mask_cols = []
    empty_mask_cols = []
    for year in cost_years:
        merged[f"Mask_grey_{year}"] = _classify_cost_points(
            merged[f"Cost_{cost_type}_{year}_orig"],
            merged["Injection Rate_merged"],
            merged["Storage Potential_merged"],
        )[2]
        mask_cols.append(f"Mask_grey_{year}")

        merged[f"Mask_empty_{year}"] = (
            merged[f"Cost_{cost_type}_{year}_orig"].isna()
            & merged["Storage Potential_merged"].gt(0)
        )
        empty_mask_cols.append(f"Mask_empty_{year}")

    # Preserve outlines for every valid grid cell, even when cost is missing/NaN.
    merged["Mask_coverage_any"] = merged["X"].notna() & merged["Y"].notna()

    if mask_cols:
        merged["Mask_outline_any"] = (
            merged[mask_cols].any(axis=1)
            | merged[empty_mask_cols].any(axis=1)
            | merged["Mask_coverage_any"]
        )
    else:
        merged["Mask_outline_any"] = merged["Mask_coverage_any"]

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

    # Optional companion coverage (e.g. DSA XY when plotting EOR-only) so empty
    # EOR regions still keep the same area boundaries as DSA / combined maps.
    coverage_paths = outline_coverage_paths
    if coverage_paths is None:
        coverage_paths = (excel1_path, excel2_path)
    elif isinstance(coverage_paths, (str, Path)):
        coverage_paths = (coverage_paths,)

    coverage_frames = []
    for path in coverage_paths:
        cov = pd.read_excel(path, usecols=["X", "Y"])
        cov = cov.dropna(subset=["X", "Y"])
        cov_box = (
            ((cov["X"] >= 108) & (cov["X"] <= 118) & (cov["Y"] >= 15) & (cov["Y"] <= 24))
            | ((cov["X"] >= 120) & (cov["X"] <= 129) & (cov["Y"] >= 25) & (cov["Y"] <= 34))
        )
        coverage_frames.append(cov.loc[cov_box, ["X", "Y"]])

    outline_xy = pd.concat(
        [merged.loc[merged["Mask_outline_any"], ["X", "Y"]], *coverage_frames],
        ignore_index=True,
    ).drop_duplicates()

    gdf_outline_any = gpd.GeoDataFrame(
        outline_xy,
        geometry=gpd.points_from_xy(outline_xy["X"], outline_xy["Y"]),
        crs="EPSG:4326",
    ).to_crs(epsg=3857)

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
        col_usd = f"Cost_{cost_type}_{year}_USD"
        col_orig = f"Cost_{cost_type}_{year}_orig"

        # Define point masks
        mask_pos, mask_red, mask_grey = _classify_cost_points(
            gdf_points[col_orig],
            gdf_points["Injection Rate_merged"],
            gdf_points["Storage Potential_merged"],
        )

        # Plot boundaries
        china_outline.boundary.plot(ax=axes[i], linewidth=1.5, color="black")
        if not china_provinces.empty:
            china_provinces.plot(ax=axes[i], facecolor="none", edgecolor="gray", linewidth=1)

        # Plot positive cost points with colormap
        if not gdf_points[mask_pos].empty:
            vals = gdf_points.loc[mask_pos, col_usd]
            # choose per-row vmax when colorbar==2, otherwise use global
            if colorbar == 2:
                # top row indices 0-2 use row1_vmax, bottom row 3-5 use row2_vmax
                vmax_use = row1_vmax if i < 3 else row2_vmax
            else:
                vmax_use = global_vmax
            plot_vals = vals.clip(upper=vmax_use)
            gdf_plot = gdf_points[mask_pos].copy()
            gdf_plot[f"{col_usd}_plot"] = plot_vals
            
            gdf_plot.plot(
                ax=axes[i],
                column=f"{col_usd}_plot",
                cmap=COST_CMAP,
                markersize=15,
                legend=False,
                vmin=global_vmin,
                vmax=vmax_use
            )

        # Plot exhausted/non-positive points (red category)
        if not gdf_points[mask_red].empty:
            gdf_points[mask_red].plot(ax=axes[i], color=COLOR_EXHAUSTED, markersize=15)

        # Plot special zero points (undetermined)
        if not gdf_points[mask_grey].empty:
            gdf_points[mask_grey].plot(ax=axes[i], color=COLOR_UNDETERMINED, markersize=15)

        # Plot the outline-only polygon for the two allowed regions, if any.
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

        # Format subplot
        bounds = china_outline.total_bounds
        axes[i].set_xlim(bounds[0], bounds[2])
        axes[i].set_ylim(bounds[1], bounds[3])
        axes[i].set_title(get_subplot_label(year), fontsize=mf["title"], fontweight="bold")
        axes[i].set_aspect('equal')
        axes[i].set_xticks([])
        axes[i].set_yticks([])

    # Add colorbars and legend
    # Tighten left margin so the map grid doesn't leave a large blank strip.
    plt.subplots_adjust(left=0.04, right=0.88, hspace=0.28, wspace=0.1)

    if colorbar == 1:
        # Unified colorbar for all panels (use 5 tick marks)
        cbar_ax = fig.add_axes([0.91, 0.12, 0.015, 0.76])  # [left, bottom, width, height]
        sm = plt.cm.ScalarMappable(cmap=COST_CMAP, norm=plt.Normalize(vmin=global_vmin, vmax=global_vmax))
        sm._A = []
        cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical")
        # Use exactly 5 tick positions and label the last as '≥vmax'
        ticks = np.linspace(global_vmin, global_vmax, 5)
        tick_labels = [f"{t:.0f}" for t in ticks[:-1]] + [f"≥{int(global_vmax)}"]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(tick_labels)
        cbar.set_label(f"Cost ({cost_unit_label})", fontsize=mf["colorbar_label"])
        cbar.ax.tick_params(labelsize=mf["colorbar_tick"])
    else:
        # Two colorbars: one aligned with the first row, one with the second row
        # Right-side vertical bars: top for first row, bottom for second row
        sm1 = plt.cm.ScalarMappable(cmap=COST_CMAP, norm=plt.Normalize(vmin=global_vmin, vmax=row1_vmax))
        sm1._A = []
        cbar_ax1 = fig.add_axes([0.91, 0.56, 0.015, 0.34])
        cbar1 = fig.colorbar(sm1, cax=cbar_ax1, orientation="vertical")
        ticks1 = np.linspace(global_vmin, row1_vmax, 5)
        tick_labels1 = [f"{t:.0f}" for t in ticks1[:-1]] + [f"≥{int(row1_vmax)}"]
        cbar1.set_ticks(ticks1)
        cbar1.set_ticklabels(tick_labels1)
        cbar1.set_label(f"Cost ({cost_unit_label})", fontsize=mf["colorbar_label"])
        cbar1.ax.tick_params(labelsize=mf["colorbar_tick"])

        sm2 = plt.cm.ScalarMappable(cmap=COST_CMAP, norm=plt.Normalize(vmin=global_vmin, vmax=row2_vmax))
        sm2._A = []
        cbar_ax2 = fig.add_axes([0.91, 0.12, 0.015, 0.34])
        cbar2 = fig.colorbar(sm2, cax=cbar_ax2, orientation="vertical")
        ticks2 = np.linspace(global_vmin, row2_vmax, 5)
        tick_labels2 = [f"{t:.0f}" for t in ticks2[:-1]] + [f"≥{int(row2_vmax)}"]
        cbar2.set_ticks(ticks2)
        cbar2.set_ticklabels(tick_labels2)
        cbar2.set_label(f"Cost ({cost_unit_label})", fontsize=mf["colorbar_label"])
        cbar2.ax.tick_params(labelsize=mf["colorbar_tick"])
    
    # Add legend for color categories in a single horizontal row
    legend_ax = fig.add_axes([0.02, 0.01, 0.88, 0.07])
    legend_ax.axis('off')

    legend_handles = [
         Line2D([0], [0], marker='o', color='none', markerfacecolor='red', markeredgecolor='red', markersize=mf["legend_marker"],
             label='RED: Exhausted storage potential'),
         Line2D([0], [0], marker='o', color='none', markerfacecolor=COLOR_UNDETERMINED, markeredgecolor=COLOR_UNDETERMINED, markersize=mf["legend_marker"],
             label='DEEP GREEN: Undetermined injection rate')
    ]

    legend_ax.legend(
        handles=legend_handles,
        loc='center',
        ncol=2,
        frameon=False,
        fontsize=mf["legend"],
        handletextpad=0.6,
        columnspacing=2.2
    )

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[OK] Plot saved to {output_path}")
