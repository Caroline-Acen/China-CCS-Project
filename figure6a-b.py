"""
Figure 6a–b: Province-level storage vs injection — DSA and EOR.

One 1×2 figure with a shared legend so both panels have the same height.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, ScalarFormatter

from config import (
    CARBON_DSA,
    CARBON_EOR,
    chart_path,
    manuscript_font_bundle,
)
from province_style import build_province_style

# Journal request: 4500 px wide; height 3000 px keeps both panels equal.
TARGET_WIDTH_PX = 4500
TARGET_HEIGHT_PX = 3000
DPI = 300
FIG_W_IN = TARGET_WIDTH_PX / DPI
FIG_H_IN = TARGET_HEIGHT_PX / DPI


def _prepare_plot_data(carbon_dsa_path, carbon_eor_path):
    carbon_dsa = pd.read_excel(carbon_dsa_path)
    carbon_eor = pd.read_excel(carbon_eor_path)

    all_provinces = sorted(
        set(carbon_dsa["Province name"]).union(set(carbon_eor["Province name"]))
    )
    rows = []
    for province in all_provinces:
        dsa_row = carbon_dsa[carbon_dsa["Province name"] == province]
        eor_row = carbon_eor[carbon_eor["Province name"] == province]
        dsa_storage = (
            dsa_row["Storage potential-DSA (Mt)"].values[0] if len(dsa_row) else 0
        )
        dsa_injection = (
            dsa_row["Injection rate capability-DSA (Mt/a) (Average)"].values[0]
            if len(dsa_row)
            else 0
        )
        eor_storage = (
            eor_row["Storage potential-EOR (Mt)"].values[0] if len(eor_row) else 0
        )
        eor_injection = (
            eor_row["Injection rate capability-EOR (Mt/a) (Average)"].values[0]
            if len(eor_row)
            else 0
        )
        if dsa_storage > 0 or eor_storage > 0:
            rows.append(
                {
                    "province": province.strip(),
                    "dsa_storage": dsa_storage,
                    "dsa_injection": dsa_injection,
                    "eor_storage": eor_storage,
                    "eor_injection": eor_injection,
                }
            )

    plot_df = pd.DataFrame(rows)
    style = build_province_style(plot_df["province"])
    return plot_df, style


def _format_storage_axis(ax):
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, prune=None))
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.ticklabel_format(style="plain", axis="x")
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")


def _scatter_panel(ax, plot_df, style, storage_col, injection_col, panel_label: str, mf: dict):
    subset = plot_df[plot_df[storage_col] > 0].sort_values("province")

    for row in subset.itertuples():
        st = style[row.province]
        ax.scatter(
            getattr(row, storage_col),
            getattr(row, injection_col),
            c=[st["color"]],
            marker=st["marker"],
            s=160,
            alpha=0.9,
            edgecolors="black",
            linewidths=0.6,
        )

    ax.set_xlabel("Storage Capacity (Mt CO$_2$)", fontsize=mf["label"], fontweight="bold")
    ax.set_ylabel("Injection Rate (Mt CO$_2$/year)", fontsize=mf["label"], fontweight="bold")
    ax.tick_params(labelsize=mf["tick"])
    ax.grid(False)
    _format_storage_axis(ax)
    ax.text(
        0.5,
        0.98,
        panel_label,
        transform=ax.transAxes,
        va="top",
        ha="center",
        fontsize=mf["title"],
        fontweight="bold",
        zorder=10,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="none"),
    )


def _add_province_legend(fig, axes, style, legend_provinces, mf: dict):
    labels = list(legend_provinces)
    legend_fs = mf["legend"]
    marker_fs = max(14.0, round(legend_fs * 0.8, 1))
    handles = [
        Line2D(
            [0],
            [0],
            marker=style[name]["marker"],
            color="w",
            markerfacecolor=style[name]["color"],
            markeredgecolor="black",
            markeredgewidth=0.65,
            markersize=marker_fs,
            linestyle="",
        )
        for name in labels
    ]
    ncol = math.ceil(len(labels) / 4)
    axes_list = list(np.atleast_1d(axes).flat)
    left = 0.01
    right = min(0.99, axes_list[-1].get_position().x1)
    fig.legend(
        handles,
        labels,
        loc="lower left",
        bbox_to_anchor=(left, 0.02, right - left, 0.16),
        mode="expand",
        ncol=ncol,
        prop={"size": legend_fs, "weight": "normal"},
        frameon=True,
        fancybox=False,
        handletextpad=0.6,
        columnspacing=0.85,
        labelspacing=0.45,
        borderpad=0.35,
        handlelength=1.7,
        borderaxespad=0.0,
    )


def plot_figure6ab(carbon_dsa_path, carbon_eor_path, output_path: str) -> None:
    plot_df, style = _prepare_plot_data(carbon_dsa_path, carbon_eor_path)
    legend_provinces = sorted(plot_df["province"].tolist())

    fig, axes = plt.subplots(1, 2, figsize=(FIG_W_IN, FIG_H_IN), dpi=DPI)
    mf = manuscript_font_bundle(fig.get_figwidth(), target_pt=9.0)
    _scatter_panel(axes[0], plot_df, style, "dsa_storage", "dsa_injection", "(a) DSA", mf)
    _scatter_panel(axes[1], plot_df, style, "eor_storage", "eor_injection", "(b) EOR", mf)

    fig.subplots_adjust(left=0.07, right=0.99, top=0.97, bottom=0.34, wspace=0.28)
    _add_province_legend(fig, axes, style, legend_provinces, mf)
    fig.savefig(output_path, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"[OK] Figure 6a–b saved to {output_path}")


if __name__ == "__main__":
    plot_figure6ab(
        carbon_dsa_path=str(CARBON_DSA),
        carbon_eor_path=str(CARBON_EOR),
        output_path=chart_path("figure6a-b.png"),
    )
