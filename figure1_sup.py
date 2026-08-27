# """
# Figure 1 Supplementary: Storage Potential and Industries Map (Original)
# China map showing all industry categories overlaid with storage potential locations.
# Similar to figure4a but used for supplementary materials.
# """

# import pandas as pd
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import matplotlib.cm as cm
# from matplotlib.lines import Line2D
# from matplotlib.patches import Patch
# from shapely.ops import unary_union
# from textwrap import wrap


# def plot_storage_potential_and_industries(
#     dsa_pipeline, 
#     eor_pipeline,
#     industries_csv, 
#     admin0_path, 
#     admin1_path,
#     output_path="storage_industry_map_sup.png"
# ):
#     """
#     Plot storage potential and industry locations with 250km buffer overlap.

#     Parameters:
#         dsa_pipeline: Path to DSA Pipeline Excel file
#         eor_pipeline: Path to EOR Pipeline Excel file
#         industries_csv: Path to industries CSV
#         admin0_path: Path to Natural Earth admin0 shapefile
#         admin1_path: Path to Natural Earth admin1 shapefile
#         output_path: Output image path
#     """

#     print("Loading data...")

#     # Load Storage Potential
#     dsa_df = pd.read_excel(dsa_pipeline)
#     eor_df = pd.read_excel(eor_pipeline)

#     dsa_gdf = gpd.GeoDataFrame(
#         dsa_df,
#         geometry=gpd.points_from_xy(dsa_df["X"], dsa_df["Y"]),
#         crs="EPSG:4326"
#     )
#     eor_gdf = gpd.GeoDataFrame(
#         eor_df,
#         geometry=gpd.points_from_xy(eor_df["X"], eor_df["Y"]),
#         crs="EPSG:4326"
#     )
    
#     # Load Industries
#     ind_df = pd.read_csv(industries_csv, low_memory=False)

#     # Industry category markers
#     categories = {
#         "TC19": "o",   # Transport
#         "ME19": "s",   # Manufacturing
#         "WF19": "^",   # Waste
#         "AF19": "D",   # Agriculture
#         "OM19": "P",   # Other industries
#         "MC19": "X",   # Mining
#         "PM19": "v",   # Paper
#         "PP19": "*",   # Power
#     }

#     # Filter rows where all categories are zero
#     ind_df = ind_df[ind_df[list(categories.keys())].sum(axis=1) > 0]
#     print(f"Filtered industries: {len(ind_df)} rows kept")

#     # Coordinate handling
#     x_min = ind_df["X"].min()
#     target_crs = "ESRI:102012"
#     already_projected = (abs(x_min) > 1000)

#     if already_projected:
#         base_gdf = gpd.GeoDataFrame(
#             ind_df,
#             geometry=gpd.points_from_xy(ind_df["X"], ind_df["Y"]),
#             crs=target_crs
#         )
#     else:
#         base_gdf = gpd.GeoDataFrame(
#             ind_df,
#             geometry=gpd.points_from_xy(ind_df["X"], ind_df["Y"]),
#             crs="EPSG:4326"
#         )

#     # Load shapefiles
#     admin0 = gpd.read_file(admin0_path)
#     admin1 = gpd.read_file(admin1_path)
    
#     china = admin0[admin0["ADMIN"] == "China"]
#     if "admin" in admin1.columns:
#         provinces = admin1[admin1["admin"] == "China"]
#     elif "ADM0NAME" in admin1.columns:
#         provinces = admin1[admin1["ADM0NAME"] == "China"]
#     else:
#         provinces = gpd.GeoDataFrame(columns=admin1.columns, geometry=admin1.geometry.name)
    
#     # Reproject all data
#     china = china.to_crs(target_crs)
#     provinces = provinces.to_crs(target_crs)
#     dsa_gdf = dsa_gdf.to_crs(target_crs)
#     eor_gdf = eor_gdf.to_crs(target_crs)
#     base_gdf = base_gdf.to_crs(target_crs)

#     # Compute 250km buffers
#     print("Computing 250km buffers...")
#     buffer_250km = base_gdf.buffer(250_000)
#     buffer_union = unary_union(buffer_250km)

#     overlap_polygons = []
#     for pt in pd.concat([dsa_gdf, eor_gdf]).geometry:
#         if buffer_union.contains(pt):
#             overlap_polygons.append(buffer_union)
    
#     overlap_gdf = gpd.GeoDataFrame(
#         geometry=[unary_union(overlap_polygons)],
#         crs=target_crs
#     ) if overlap_polygons else gpd.GeoDataFrame(geometry=[], crs=target_crs)

#     # Create figure
#     print("Creating plot...")
#     fig, ax = plt.subplots(1, 1, figsize=(10, 10))

#     china.boundary.plot(ax=ax, color="black", linewidth=1)
#     provinces.boundary.plot(ax=ax, color="gray", linewidth=0.5, alpha=0.7)

#     # Color scale
#     max_val = base_gdf[list(categories.keys())].max().max()
#     bins = [0, 1, 2, 3, 4, 5, max_val]
#     labels = ["1", "2", "3", "4", "5", f">{5}"]
#     cmap = plt.get_cmap("YlOrRd", len(labels))
#     norm = mcolors.BoundaryNorm(boundaries=bins, ncolors=len(labels), clip=False)

#     # Plot industries by category
#     for cat, marker in categories.items():
#         subset = base_gdf[base_gdf[cat] > 0]
#         if subset.empty:
#             continue
#         subset.plot(
#             ax=ax,
#             markersize=40,
#             marker=marker,
#             c=subset[cat],
#             cmap=cmap,
#             norm=norm,
#             alpha=0.7
#         )

#     # Plot storage potentials
#     dsa_gdf.plot(ax=ax, color="green", markersize=25, alpha=0.7)
#     eor_gdf.plot(ax=ax, color="gray", markersize=25, alpha=0.7)

#     # Colorbar
#     sm = cm.ScalarMappable(cmap=cmap, norm=norm)
#     sm.set_array([])
#     cbar = fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.03, pad=0.02)
#     cbar.set_ticks([(bins[i] + bins[i+1]) / 2 for i in range(len(labels))])
#     cbar.set_ticklabels(labels)
#     cbar.set_label("Number of Industries", fontsize=10)

#     # Overlap region
#     if not overlap_gdf.empty:
#         overlap_gdf.plot(ax=ax, facecolor="none", edgecolor="blue",
#                          hatch="///", linewidth=0, alpha=0.5)

#     # Legend
#     legend_elements = []
#     for cat, marker in categories.items():
#         legend_elements.append(Line2D([0], [0], marker=marker, color="w",
#                                       markerfacecolor="orange", markersize=8,
#                                       label=cat.replace("19", "")))
#     legend_elements.extend([
#         Patch(facecolor="green", edgecolor="black", label="DSA storage potential"),
#         Patch(facecolor="gray", edgecolor="black", label="EOR storage potential"),
#         Patch(facecolor="none", edgecolor="blue", hatch="///",
#               label="Industries and storage potential < 250km"),
#     ])
#     labels_wrapped = ["\n".join(wrap(elem.get_label(), 30)) for elem in legend_elements]
#     ax.legend(handles=legend_elements, labels=labels_wrapped,
#               loc="lower left", fontsize=9, frameon=True)

#     ax.axis("off")
    
#     plt.savefig(output_path, dpi=300, bbox_inches="tight")
#     plt.close()
#     print(f"[OK] Map saved to {output_path}")


# # Run function
# plot_storage_potential_and_industries(
#     dsa_pipeline="./data/DSA-Pipeline.xlsx",
#     eor_pipeline="./data/EOR_Pipeline.xlsx",
#     industries_csv="./data/industries.csv",
#     admin0_path="./data/natural_earth/admin0",
#     admin1_path="./data/natural_earth/admin1",
#     output_path="storage_industry_map_sup.png"
# )





"""
Figure 1: Injection Maps
Plots China maps showing Storage Potential and Injection Rate from 4 Excel files.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from config import DSA_PIPELINE, DSA_TANK, EOR_PIPELINE, EOR_TANK, NE_ADMIN0, NE_ADMIN1, chart_path, manuscript_font_bundle
from cost_maps import COST_CMAP

# Same sequential scheme as cost maps: deep blue → cyan → yellow → brown
INJECTION_CMAP = COST_CMAP


def plot_injection_maps(
    excel1, excel2, excel3, excel4,
    admin0_path, admin1_path,
    output_path="china_injection_maps.png",
    colorbar=2
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
        - COST_CMAP: Positive values (blue → cyan → yellow → brown, same as cost maps)

    Parameter `colorbar`:
        - 2 (default): two vertical colorbars, one for the first row and one for the second row
        - 4: one vertical colorbar for each subplot
    """

    if colorbar not in [2, 4]:
        raise ValueError("`colorbar` must be either 2 or 4.")

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
    titles = ["(a)", "(b)", "(c)", "(d)"]

    # Layout: for colorbar=4 use map | cbar | spacer | map | cbar so labels never hit the next map
    if colorbar == 4:
        fig = plt.figure(figsize=(20, 14), constrained_layout=False)
        gs = fig.add_gridspec(
            2, 5,
            # Spacer keeps ~20% of prior gap; maps reclaim the other 80%
            width_ratios=[1.128, 0.05, 0.104, 1.128, 0.05],
            height_ratios=[1.0, 1.0],
            left=0.03, right=0.98, top=0.95, bottom=0.05,
            wspace=0.18, hspace=0.18,
        )
        axes = np.array([
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 3]),
            fig.add_subplot(gs[1, 0]),
            fig.add_subplot(gs[1, 3]),
        ])
        caxes = [
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[0, 4]),
            fig.add_subplot(gs[1, 1]),
            fig.add_subplot(gs[1, 4]),
        ]
    else:
        fig, axes = plt.subplots(2, 2, figsize=(18, 14), constrained_layout=False)
        axes = axes.ravel()
        caxes = None

    mf = manuscript_font_bundle(fig.get_figwidth())
    fs_title = mf["title"] * 1.15
    fs_cbar = mf["colorbar_label"] * 1.15
    fs_tick = mf["colorbar_tick"] * 1.15

    # Compute a global vmax across all subplots (use per-subplot percentile then combine)
    all_vals = []
    for i, col in enumerate(columns):
        vals = gdfs[i][gdfs[i][col] > 0][col]
        if not vals.empty:
            perc = subplot_percentiles.get(i, 0.95)
            vmax_i = vals.quantile(perc)
            clipped = vals.clip(upper=vmax_i)
            all_vals.extend(clipped.dropna().tolist())

    if all_vals:
        global_vmin = 0
        global_vmax = float(np.percentile(all_vals, 95))
    else:
        global_vmin, global_vmax = 0, 1

    # If colorbar==2, compute separate vmax for top row (indices 0-1) and bottom row (2-3)
    if colorbar == 2:
        row1_vals = []
        row2_vals = []
        for i, col in enumerate(columns):
            vals = gdfs[i][gdfs[i][col] > 0][col]
            if vals.empty:
                continue
            perc = subplot_percentiles.get(i, 0.95)
            vmax_i = vals.quantile(perc)
            clipped = vals.clip(upper=vmax_i)
            if i < 2:
                row1_vals.extend(clipped.dropna().tolist())
            else:
                row2_vals.extend(clipped.dropna().tolist())

        row1_vmax = float(np.percentile(row1_vals, 95)) if row1_vals else global_vmax
        row2_vmax = float(np.percentile(row2_vals, 95)) if row2_vals else global_vmax
    else:
        row1_vmax = row2_vmax = global_vmax

    # Per-subplot vmax values for colorbar==4
    subplot_vmax = []
    for i, col in enumerate(columns):
        vals = gdfs[i][gdfs[i][col] > 0][col]
        if vals.empty:
            subplot_vmax.append(global_vmax)
        else:
            perc = subplot_percentiles.get(i, 0.95)
            subplot_vmax.append(float(vals.quantile(perc)))

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
            vals = gdf[mask_pos][col]

            if colorbar == 2:
                vmax_use = row1_vmax if i < 2 else row2_vmax
            else:
                vmax_use = subplot_vmax[i]
            
            plot_vals = vals.clip(upper=vmax_use)
            gdf_plot = gdf[mask_pos].copy()
            gdf_plot[f"{col}_plot"] = plot_vals
            
            gdf_plot.plot(
                ax=ax,
                column=f"{col}_plot",
                cmap=INJECTION_CMAP,
                markersize=15,
                legend=False,
                vmin=global_vmin,
                vmax=vmax_use
            )

        # Plot negative/zero points (red exhausted category)
        if not gdf[mask_red].empty:
            gdf[mask_red].plot(ax=ax, color="red", markersize=15)

        # Plot special zero points (grey)
        if not gdf[mask_grey].empty:
            gdf[mask_grey].plot(ax=ax, color="grey", markersize=15)

        # Format subplot
        bounds = china_outline.total_bounds
        ax.set_xlim(bounds[0], bounds[2])
        ax.set_ylim(bounds[1], bounds[3])
        ax.set_title(title, fontsize=fs_title, fontweight="bold")
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    # Colorbar layout
    def _style_cbar(cbar, label: str) -> None:
        cbar.locator = MaxNLocator(nbins=6)
        cbar.update_ticks()
        # labelpad keeps title clear of tick numbers; spacer col absorbs the rest
        cbar.set_label(label, fontsize=fs_cbar, labelpad=8)
        cbar.ax.tick_params(labelsize=fs_tick, pad=2)
        cbar.outline.set_linewidth(0.6)

    if colorbar == 2:
        left_margin = 0.04
        right_margin = 0.86
        fig.subplots_adjust(
            left=left_margin, right=right_margin, top=0.96, bottom=0.05,
            hspace=0.14, wspace=0.10,
        )
        fig.canvas.draw()
        top_axes = [axes[0], axes[1]]
        bot_axes = [axes[2], axes[3]]
        y0_top = min(ax.get_position().y0 for ax in top_axes)
        y1_top = max(ax.get_position().y1 for ax in top_axes)
        y0_bot = min(ax.get_position().y0 for ax in bot_axes)
        y1_bot = max(ax.get_position().y1 for ax in bot_axes)
        cbar_x = max(ax.get_position().x1 for ax in axes) + 0.02
        cbar_w = 0.018

        sm_top = plt.cm.ScalarMappable(
            cmap=INJECTION_CMAP, norm=plt.Normalize(vmin=global_vmin, vmax=row1_vmax),
        )
        sm_top.set_array([])
        cbar_top = fig.colorbar(
            sm_top,
            cax=fig.add_axes([cbar_x, y0_top, cbar_w, y1_top - y0_top]),
            orientation="vertical",
        )
        _style_cbar(cbar_top, "Storage Potential (Mt)")

        sm_bot = plt.cm.ScalarMappable(
            cmap=INJECTION_CMAP, norm=plt.Normalize(vmin=global_vmin, vmax=row2_vmax),
        )
        sm_bot.set_array([])
        cbar_bot = fig.colorbar(
            sm_bot,
            cax=fig.add_axes([cbar_x, y0_bot, cbar_w, y1_bot - y0_bot]),
            orientation="vertical",
        )
        _style_cbar(cbar_bot, "Injection Rate (Mt/year)")
    else:
        # Dedicated colorbar axes already allocated in the gridspec (with spacer column)
        for i, (cax, col) in enumerate(zip(caxes, columns)):
            label = "Storage Potential (Mt)" if col == "Storage Potential" else "Injection Rate (Mt/year)"
            sm = plt.cm.ScalarMappable(
                cmap=INJECTION_CMAP,
                norm=plt.Normalize(vmin=global_vmin, vmax=subplot_vmax[i]),
            )
            sm.set_array([])
            cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
            _style_cbar(cbar, label)

    plt.savefig(output_path, dpi=600, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    print(f"[OK] Plot saved to {output_path}")


# Run function
plot_injection_maps(
    excel1=str(DSA_PIPELINE),
    excel2=str(EOR_PIPELINE),
    excel3=str(DSA_TANK),
    excel4=str(EOR_TANK),
    admin0_path=str(NE_ADMIN0),
    admin1_path=str(NE_ADMIN1),
    output_path=chart_path("injection_maps.png"),
    colorbar=4,
)


