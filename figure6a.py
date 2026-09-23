"""
Figure 6a–b: Province-level storage vs injection — DSA and EOR.

Combined 1×2 subplots plus optional separate panels for readability.
"""

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


def _prepare_plot_data(carbon_dsa_path, carbon_eor_path):
    carbon_dsa = pd.read_excel(carbon_dsa_path)
    carbon_eor = pd.read_excel(carbon_eor_path)

    all_provinces = sorted(set(carbon_dsa["Province name"]).union(set(carbon_eor["Province name"])))
    rows = []
    for province in all_provinces:
        dsa_row = carbon_dsa[carbon_dsa["Province name"] == province]
        eor_row = carbon_eor[carbon_eor["Province name"] == province]
        dsa_storage = dsa_row["Storage potential-DSA (Mt)"].values[0] if len(dsa_row) else 0
        dsa_injection = dsa_row["Injection rate capability-DSA (Mt/a) (Average)"].values[0] if len(dsa_row) else 0
        eor_storage = eor_row["Storage potential-EOR (Mt)"].values[0] if len(eor_row) else 0
        eor_injection = eor_row["Injection rate capability-EOR (Mt/a) (Average)"].values[0] if len(eor_row) else 0
        if dsa_storage > 0 or eor_storage > 0:
            rows.append({
                "province": province.strip(),
                "dsa_storage": dsa_storage,
                "dsa_injection": dsa_injection,
                "eor_storage": eor_storage,
                "eor_injection": eor_injection,
            })

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
        0.5, 0.98, panel_label,
        transform=ax.transAxes, va="top", ha="center",
        fontsize=mf["title"], fontweight="bold",
        zorder=10,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="none"),
    )


def _add_province_legend(fig, axes, style, legend_provinces, mf: dict):
    labels = list(legend_provinces)
    legend_fs = mf["legend"]
    marker_fs = max(14.0, round(legend_fs * 0.8, 1))
    handles = [
        Line2D(
            [0], [0],
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
    # Reach nearly to the figure's left edge (past the y-label margin).
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


def plot_storage_comparison_combined(
    plot_df,
    style,
    output_path,
    *,
    legend_provinces=None,
):
    # Same canvas as the pre-/4 layout; gap under x-labels is half of that version
    fig, axes = plt.subplots(1, 2, figsize=(18, 11), dpi=600)
    mf = manuscript_font_bundle(fig.get_figwidth(), target_pt=9.0)
    _scatter_panel(axes[0], plot_df, style, "dsa_storage", "dsa_injection", "(a) DSA", mf)
    _scatter_panel(axes[1], plot_df, style, "eor_storage", "eor_injection", "(b) EOR", mf)

    if legend_provinces:
        # Pre-/4 used bottom=0.46 (gap ≈ 0.25); half that gap → bottom ≈ 0.334
        fig.subplots_adjust(left=0.06, right=0.99, top=0.97, bottom=0.334, wspace=0.28)
        _add_province_legend(fig, axes, style, legend_provinces, mf)
        fig.savefig(output_path, dpi=600, facecolor="white")
    else:
        fig.tight_layout()
        fig.savefig(output_path, dpi=600, bbox_inches="tight", pad_inches=0.15)

    plt.close(fig)
    print(f"[OK] Combined plot saved to {output_path}")


def _plot_single_panel(
    plot_df,
    style,
    storage_col,
    injection_col,
    output_path,
    *,
    panel_label: str = "",
    show_legend=False,
    legend_provinces=None,
):
    fig, ax = plt.subplots(figsize=(11, 9), dpi=600)
    mf = manuscript_font_bundle(fig.get_figwidth(), target_pt=9.0)
    _scatter_panel(ax, plot_df, style, storage_col, injection_col, panel_label, mf)

    if show_legend and legend_provinces:
        fig.subplots_adjust(left=0.10, right=0.98, top=0.96, bottom=0.334)
        _add_province_legend(fig, ax, style, legend_provinces, mf)
        fig.savefig(output_path, dpi=600, facecolor="white")
    else:
        fig.tight_layout()
        fig.savefig(output_path, dpi=600, bbox_inches="tight", pad_inches=0.15)

    plt.close(fig)
    print(f"[OK] Plot saved to {output_path}")


def plot_storage_comparison(carbon_dsa_path, carbon_eor_path, output_dir=None):
    plot_df, style = _prepare_plot_data(carbon_dsa_path, carbon_eor_path)
    legend_provinces = sorted(plot_df["province"].tolist())

    combined_path = (
        chart_path("storage_comparison.png")
        if output_dir is None
        else f"{output_dir}/storage_comparison.png"
    )
    dsa_path = (
        chart_path("storage_comparison_dsa.png")
        if output_dir is None
        else f"{output_dir}/storage_comparison_dsa.png"
    )
    eor_path = (
        chart_path("storage_comparison_eor.png")
        if output_dir is None
        else f"{output_dir}/storage_comparison_eor.png"
    )

    plot_storage_comparison_combined(
        plot_df, style, combined_path, legend_provinces=legend_provinces,
    )
    _plot_single_panel(
        plot_df, style, "dsa_storage", "dsa_injection", dsa_path,
        panel_label="(a) DSA",
    )
    _plot_single_panel(
        plot_df, style, "eor_storage", "eor_injection", eor_path,
        panel_label="(b) EOR",
        show_legend=True,
        legend_provinces=legend_provinces,
    )


if __name__ == "__main__":
    plot_storage_comparison(
        carbon_dsa_path=str(CARBON_DSA),
        carbon_eor_path=str(CARBON_EOR),
    )
