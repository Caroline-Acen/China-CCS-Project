"""
Figure 5c (submitted): Province-level storage vs injection — DSA and EOR.

Combined 1×2 subplots plus optional separate panels for readability.
"""

import math

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


def _add_province_legend(fig, style, legend_provinces, mf: dict):
    labels = list(legend_provinces)
    handles = [
        Line2D(
            [0], [0],
            marker=style[name]["marker"],
            color="w",
            markerfacecolor=style[name]["color"],
            markeredgecolor="black",
            markeredgewidth=0.6,
            markersize=mf["legend_marker"],
            linestyle="",
        )
        for name in labels
    ]
    n_items = len(labels)
    ncol = math.ceil(n_items / 3)
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.02, 0.01, 0.96, 0.14),
        mode="expand",
        ncol=ncol,
        fontsize=mf["legend"],
        frameon=True,
        handletextpad=0.35,
        columnspacing=1.2,
        labelspacing=0.5,
    )


def plot_storage_comparison_combined(
    plot_df,
    style,
    output_path,
    *,
    legend_provinces=None,
):
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=600)
    mf = manuscript_font_bundle(fig.get_figwidth())
    _scatter_panel(axes[0], plot_df, style, "dsa_storage", "dsa_injection", "(a) DSA", mf)
    _scatter_panel(axes[1], plot_df, style, "eor_storage", "eor_injection", "(b) EOR", mf)

    if legend_provinces:
        _add_province_legend(fig, style, legend_provinces, mf)
        fig.tight_layout(rect=[0, 0.16, 1, 1])
    else:
        fig.tight_layout()

    fig.savefig(output_path, dpi=600, bbox_inches="tight")
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
    mf = manuscript_font_bundle(fig.get_figwidth())
    _scatter_panel(ax, plot_df, style, storage_col, injection_col, panel_label, mf)

    if show_legend and legend_provinces:
        _add_province_legend(fig, style, legend_provinces, mf)
        fig.tight_layout(rect=[0, 0.18, 1, 1])
    else:
        fig.tight_layout()

    fig.savefig(output_path, dpi=600, bbox_inches="tight")
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
