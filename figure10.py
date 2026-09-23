"""
Figure 10: Responsibility vs capability bubble charts (DSA and EOR).
Reviewer fix: true 1:1 reference line on log-log axes.
"""

import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext, LogLocator

from config import (
    CARBON_DSA,
    CARBON_EOR,
    CHARTS_DIR,
    ensure_charts_dir,
    manuscript_font_bundle,
)
from province_style import province_style_from_carbon


def _apply_log_axes(ax):
    """Log-log axes with mathtext power labels (10^n) on both axes."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10))
        axis.set_major_formatter(LogFormatterMathtext(base=10))


def _add_one_to_one_line(ax, x_vals, y_vals):
    positive = (x_vals > 0) & (y_vals > 0)
    if positive.sum() < 2:
        return
    xmin = x_vals[positive].min()
    xmax = x_vals[positive].max()
    ymin = y_vals[positive].min()
    ymax = y_vals[positive].max()
    low = min(xmin, ymin)
    high = max(xmax, ymax)
    line = np.logspace(np.log10(low), np.log10(high), 100)
    ax.plot(line, line, color="red", alpha=0.45, linewidth=2, zorder=1)


def generate_bubble_charts(dsa_file, eor_file, output_dir="."):
    dsa_df = pd.read_excel(dsa_file)
    eor_df = pd.read_excel(eor_file)
    dsa_df["Province name"] = dsa_df["Province name"].str.strip()
    eor_df["Province name"] = eor_df["Province name"].str.strip()
    style = province_style_from_carbon(dsa_file, eor_file)
    years = ["2025", "2030", "2035", "2040", "2045", "2050"]

    def create_bubble_plot(df, storage_type):
        fig, axes = plt.subplots(2, 3, figsize=(20, 13))
        mf = manuscript_font_bundle(fig.get_figwidth())
        axes = axes.flatten()
        injection_col = f"Injection rate capability-{storage_type} (Mt/a) (Average)"
        injection_rates = df[injection_col].values
        rate_range = injection_rates.max() - injection_rates.min()
        bubble_sizes = (
            (injection_rates - injection_rates.min()) / rate_range * 500 + 50
            if rate_range > 0 else np.full_like(injection_rates, 100, dtype=float)
        )

        for i, year in enumerate(years):
            ax = axes[i]
            emission_col = f"CO2_{year}"
            storage_col = f"P_{year}"
            mask = (df[storage_col] > 0) & (df[emission_col] > 0)
            plot_data = df[mask].copy()

            x_vals = plot_data[emission_col].to_numpy()
            y_vals = plot_data[storage_col].to_numpy()

            for province in plot_data["Province name"]:
                row = plot_data[plot_data["Province name"] == province].iloc[0]
                idx = row.name
                st = style.get(province, {"color": "#333333", "marker": "o"})
                ax.scatter(
                    row[emission_col],
                    row[storage_col],
                    s=bubble_sizes[idx],
                    c=[st["color"]],
                    marker=st["marker"],
                    alpha=0.85,
                    edgecolors="black",
                    linewidths=0.4,
                    label=province if i == 0 else "",
                    zorder=3,
                )

            if len(plot_data) > 1:
                _apply_log_axes(ax)
                _add_one_to_one_line(ax, x_vals, y_vals)

            ax.set_xlabel("CO$_2$ emissions (Mt)", fontsize=mf["label"], fontweight="bold")
            ax.set_ylabel("Storage potential (Mt)", fontsize=mf["label"], fontweight="bold")
            ax.set_title(year, fontsize=mf["title"], fontweight="bold")
            ax.tick_params(axis="both", labelsize=mf["tick"])
            ax.grid(False)

        handles, labels = axes[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        has_parity = True

        legend_labels = sorted(by_label.keys())
        uniform_handles = [
            Line2D(
                [0], [0],
                marker=style.get(name, {"marker": "o"})["marker"],
                color="w",
                markerfacecolor=style.get(name, {"color": "#333333"})["color"],
                markeredgecolor="black",
                markeredgewidth=0.4,
                markersize=mf["legend_marker"],
                linestyle="",
                alpha=0.85,
            )
            for name in legend_labels
        ]
        if has_parity:
            uniform_handles.append(
                Line2D([0], [0], color="red", alpha=0.45, linewidth=2)
            )
            legend_labels.append("1:1 parity")

        n_items = len(legend_labels)
        # Spread provinces across ~3 rows so markers do not overlap
        legend_rows = 3
        ncol = math.ceil(n_items / legend_rows)

        fig.legend(
            uniform_handles,
            legend_labels,
            loc="lower center",
            bbox_to_anchor=(0.01, 0.002, 0.98, 0.16),
            mode="expand",
            ncol=ncol,
            fontsize=mf["legend"],
            frameon=True,
            fancybox=False,
            handletextpad=0.35,
            columnspacing=1.4,
            labelspacing=0.55,
            borderaxespad=0.2,
            markerscale=1.0,
        )
        fig.tight_layout(rect=[0, 0.18, 1, 1])
        fig.subplots_adjust(bottom=0.18, hspace=0.32, left=0.06, right=0.98)
        return fig

    print("Generating DSA bubble chart...")
    dsa_fig = create_bubble_plot(dsa_df, "DSA")
    dsa_path = f"{output_dir}/Carbon_DSA_bubble_chart.png"
    dsa_fig.savefig(dsa_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(dsa_fig)
    print(f"[OK] DSA chart saved to: {dsa_path}")

    print("Generating EOR bubble chart...")
    eor_fig = create_bubble_plot(eor_df, "EOR")
    eor_path = f"{output_dir}/Carbon_EOR_bubble_chart.png"
    eor_fig.savefig(eor_path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(eor_fig)
    print(f"[OK] EOR chart saved to: {eor_path}")


if __name__ == "__main__":
    ensure_charts_dir()
    generate_bubble_charts(str(CARBON_DSA), str(CARBON_EOR), output_dir=str(CHARTS_DIR))
