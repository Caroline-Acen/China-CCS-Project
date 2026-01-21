"""
Carbon Contributions Module
===========================
Functions for creating carbon contribution bar charts by province and year.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def create_carbon_contribution_barcharts(
    file_path,
    output_path="Carbon_Contributions_All_Provinces.png",
    year_cols=('C_2025', 'C_2030', 'C_2035', 'C_2040', 'C_2045', 'C_2050'),
    year_labels=('2025', '2030', '2035', '2040', '2045', '2050'),
    per_province_height=0.34,
    min_subplot_height=6.0,
    left_margin=0.28,
    dpi=600
):
    """
    Create 6 horizontal bar charts (2x3 grid) showing Carbon Contributions by province.

    Parameters:
        file_path: Path to Excel file with province carbon data
        output_path: Output image path
        year_cols: Column names for each year's carbon contribution
        year_labels: Display labels for years
        per_province_height: Height per province (inches)
        min_subplot_height: Minimum subplot height (inches)
        left_margin: Space for province name labels
        dpi: Output resolution

    Returns:
        Dictionary with analysis results
    """

    # Load data
    df = pd.read_excel(file_path)
    provinces = df['Province name'].astype(str).tolist()
    n_prov = len(provinces)

    # Create consistent color mapping for provinces
    if n_prov <= 20:
        base_colors = plt.cm.tab20(np.linspace(0, 1, 20))
        colors = base_colors[:n_prov]
    else:
        colors = plt.cm.hsv(np.linspace(0, 1, n_prov, endpoint=False))
    province_to_color = {p: colors[i] for i, p in enumerate(provinces)}

    # Calculate figure dimensions
    subplot_h = max(min_subplot_height, per_province_height * n_prov)
    nrows, ncols = 2, 3
    fig_w = 22
    fig_h = subplot_h * nrows
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), dpi=dpi)
    axes = axes.flatten()

    # Font sizes
    title_fs = 18
    label_fs = 16
    tick_fs = 14
    barlabel_fs = 12

    # Plot each year
    for i, (year_col, year_label) in enumerate(zip(year_cols, year_labels)):
        ax = axes[i]

        # Get values and sort by ascending value
        vals = df.set_index('Province name')[year_col].reindex(provinces).astype(float).tolist()
        order = np.argsort(vals)
        vals_sorted = [vals[k] for k in order]
        prov_sorted = [provinces[k] for k in order]
        cols_sorted = [province_to_color[p] for p in prov_sorted]

        # Plot horizontal bars
        y_pos = np.arange(n_prov)
        bars = ax.barh(y_pos, vals_sorted, color=cols_sorted, alpha=0.9, edgecolor='white', linewidth=0.6)

        # Add value labels
        vmax = max(vals_sorted) if vals_sorted else 0.0
        x_shift = (0.01 * vmax) if vmax > 0 else 0.01
        for bar, v in zip(bars, vals_sorted):
            if v == 0:
                txt = '0'
            elif v < 0.01:
                txt = f'{v:.3f}'
            elif v < 0.1:
                txt = f'{v:.2f}'
            else:
                txt = f'{v:.1f}'
            ax.text(
                bar.get_width() + x_shift,
                bar.get_y() + bar.get_height() / 2,
                txt,
                va='center', ha='left',
                fontsize=barlabel_fs, fontweight='bold'
            )

        # Format axes
        ax.set_yticks(y_pos)
        ax.set_yticklabels(prov_sorted, fontsize=tick_fs, fontweight='bold')
        ax.set_xlabel('Carbon Contributions', fontsize=label_fs, fontweight='bold')
        
        for tick in ax.get_xticklabels():
            tick.set_fontsize(tick_fs)
            tick.set_fontweight('bold')

        ax.set_title(year_label, fontsize=title_fs, fontweight='bold')
        ax.grid(False)
        
        for spine in ax.spines.values():
            spine.set_visible(False)

        ymax = vmax * 1.10 if vmax > 0 else 1.0
        ax.set_xlim(0, ymax)

    # Remove unused axes
    for k in range(len(year_cols), len(axes)):
        fig.delaxes(axes[k])

    # Adjust layout
    plt.subplots_adjust(
        left=left_margin,
        right=0.98,
        top=0.95,
        bottom=0.08,
        wspace=0.5,
        hspace=0.25
    )

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()

    # Print summary
    print("CARBON CONTRIBUTIONS ANALYSIS")
    print("=" * 50)
    print(f"Total provinces analyzed: {n_prov}")
    totals = df[list(year_cols)].sum()
    print("\nTotal Carbon Contributions by Year:")
    for yc, yl in zip(year_cols, year_labels):
        print(f"  {yl}: {totals[yc]:.2f}")
    print(f"\n[OK] Chart saved to {output_path}")

    return {
        "n_provinces": n_prov,
        "year_totals": totals,
        "output_file": output_path
    }
