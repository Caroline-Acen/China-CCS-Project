"""
Figure 11: CCS component unit costs — China, USA, and EU.

Supplementary Table 15. Stacked bars show each component's min–max range
(height = max − min), starting at capture min.

EU values are converted with 1 USD = 0.87 EUR. USA and EU transport are
already USD/ton or EUR/ton. China transport is USD t⁻¹ km⁻¹ and is scaled
by the model source–sink haul (mean or median Best_Distance_km), not the
250 km matching radius.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from config import (
    MANUSCRIPT_FIG_WIDTH_IN,
    SS_SCORES_DIR,
    chart_path,
    data_path,
    manuscript_font_bundle,
)

CM = 1 / 2.54
FIG_W_IN = 8.5 * CM
FIG_H_IN = 9.2 * CM

EUR_PER_USD = 0.87
USD_PER_EUR = 1.0 / EUR_PER_USD
CHINA_TRANSPORT_PER_KM = (0.11, 0.19)
MATCHING_RADIUS_KM = 250.0
SUMMARY_CSV = Path(data_path("sink_suitability_summary.csv"))
HAUL_CSV = Path(data_path("source_sink_model_haul_km.csv"))
COMBINED_SCORES = SS_SCORES_DIR / "sink_suitability_combined.csv"

CAPTURE = "#4C78A8"
TRANSPORT = "#F58518"
STORAGE = "#59A14F"


def _eur_to_usd(lo: float, hi: float) -> tuple[float, float]:
    return (lo * USD_PER_EUR, hi * USD_PER_EUR)


def _per_km_to_per_ton(lo: float, hi: float, haul_km: float) -> tuple[float, float]:
    return (lo * haul_km, hi * haul_km)


def component_data(china_haul_km: float) -> dict[str, dict[str, tuple[float, float]]]:
    return {
        "China": {
            "Capture": (32.02, 43.16),
            "Transportation": _per_km_to_per_ton(*CHINA_TRANSPORT_PER_KM, china_haul_km),
            "Storage": (6.96, 8.35),
        },
        "USA": {
            "Capture": (25.0, 30.0),
            "Transportation": (25.0, 35.0),
            "Storage": (25.0, 35.0),
        },
        "EU": {
            "Capture": _eur_to_usd(40.0, 90.0),
            "Transportation": _eur_to_usd(2.0, 30.0),
            "Storage": _eur_to_usd(5.0, 35.0),
        },
    }


# Default DATA uses the 250 km radius (legacy). Mean/median variants pass data=.
DATA = component_data(MATCHING_RADIUS_KM)


def yearly_haul_table(summary_path: Path | None = None) -> pd.DataFrame:
    path = summary_path or SUMMARY_CSV
    df = pd.read_csv(path)
    keep = [
        "year",
        "n_grids",
        "best_distance_km_mean",
        "best_distance_km_median",
        "best_distance_km_p85",
        "best_distance_km_max",
    ]
    return df[keep].copy()


def pooled_mean_km(yearly: pd.DataFrame) -> float:
    w = yearly["n_grids"].to_numpy(dtype=np.float64)
    m = yearly["best_distance_km_mean"].to_numpy(dtype=np.float64)
    return float(np.average(m, weights=w))


def pooled_median_km_from_combined(
    combined_path: Path | None = None,
    chunksize: int = 500_000,
) -> tuple[float, float, int]:
    """Return (median_km, mean_km, n) from Best_Distance_km in the combined SS file."""
    path = combined_path or COMBINED_SCORES
    if not path.exists():
        raise FileNotFoundError(path)

    bins = np.linspace(0.0, MATCHING_RADIUS_KM, 25_001)
    counts = np.zeros(len(bins) - 1, dtype=np.int64)
    total = 0
    total_sum = 0.0
    for chunk in pd.read_csv(path, usecols=["Best_Distance_km"], chunksize=chunksize):
        d = pd.to_numeric(chunk["Best_Distance_km"], errors="coerce").to_numpy()
        d = d[np.isfinite(d) & (d > 0)]
        if d.size == 0:
            continue
        total += int(d.size)
        total_sum += float(d.sum())
        counts += np.histogram(d, bins=bins)[0]

    if total == 0:
        raise RuntimeError(f"No Best_Distance_km values in {path}")

    cdf = np.cumsum(counts)
    median_bin = int(np.searchsorted(cdf, total / 2.0, side="left"))
    median_bin = min(median_bin, len(counts) - 1)
    median = 0.5 * (bins[median_bin] + bins[median_bin + 1])
    mean = total_sum / total
    return float(median), float(mean), total


def write_haul_csv(
    yearly: pd.DataFrame,
    overall_mean: float,
    overall_median: float,
    overall_n: int,
    output_path: Path | None = None,
) -> Path:
    out = output_path or HAUL_CSV
    rows = yearly.copy()
    rows.insert(0, "scope", "year")
    overall = pd.DataFrame(
        [
            {
                "scope": "all_years_pooled",
                "year": "",
                "n_grids": overall_n,
                "best_distance_km_mean": overall_mean,
                "best_distance_km_median": overall_median,
                "best_distance_km_p85": np.nan,
                "best_distance_km_max": yearly["best_distance_km_max"].max(),
            }
        ]
    )
    table = pd.concat([rows, overall], ignore_index=True)
    table.to_csv(out, index=False)
    return out


def plot_figure11(output_path: str, data: dict | None = None) -> None:
    data = data or DATA
    countries = ("China", "USA", "EU")
    components = (
        ("Capture", CAPTURE),
        ("Transportation", TRANSPORT),
        ("Storage", STORAGE),
    )
    bar_w = 0.45

    fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
    mf = manuscript_font_bundle(MANUSCRIPT_FIG_WIDTH_IN)
    fs = mf["body"]

    ymax = 0.0
    for i, country in enumerate(countries):
        bottom = data[country]["Capture"][0]
        for comp, color in components:
            lo, hi = data[country][comp]
            height = hi - lo
            ax.bar(
                i,
                height,
                bottom=bottom,
                width=bar_w,
                color=color,
                edgecolor="black",
                linewidth=0.5,
                zorder=3,
            )
            bottom += height
        ymax = max(ymax, bottom)

    for x in (0.5, 1.5):
        ax.axvline(x, color="0.7", linestyle="--", linewidth=0.8, zorder=1)

    ax.set_xticks(range(len(countries)))
    ax.set_xticklabels(list(countries))
    ax.set_ylabel("Cost (USD/ton)")
    ax.set_xlabel("Countries")
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(0, ymax * 1.08)
    ax.tick_params(axis="both", labelsize=fs, width=0.6, length=3)
    ax.xaxis.label.set_size(fs)
    ax.yaxis.label.set_size(fs)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontsize(fs)

    fig.legend(
        handles=[
            Patch(facecolor=CAPTURE, edgecolor="black", linewidth=0.5, label="Capture"),
            Patch(facecolor=TRANSPORT, edgecolor="black", linewidth=0.5, label="Transportation"),
            Patch(facecolor=STORAGE, edgecolor="black", linewidth=0.5, label="Storage"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=3,
        fontsize=mf["legend"],
        frameon=False,
        handlelength=1.2,
        columnspacing=1.2,
        borderaxespad=0.0,
    )

    fig.subplots_adjust(left=0.22, right=0.96, top=0.97, bottom=0.20)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved Figure 11 to {output_path}")


def _export_high_res(src: str) -> None:
    from export_high_resolution_charts import export_one, high_res_dir

    src_path = Path(src)
    print(export_one(src_path, high_res_dir() / src_path.name))


def main() -> None:
    yearly = yearly_haul_table()
    weighted_mean = pooled_mean_km(yearly)
    print("Yearly best source–sink distance (km)")
    print(yearly.to_string(index=False))
    print(f"n-weighted mean of yearly means: {weighted_mean:.3f} km")

    median_km, combined_mean_km, n = pooled_median_km_from_combined()
    print(
        f"Pooled Best_Distance_km in combined SS file: "
        f"n={n:,}; mean={combined_mean_km:.3f} km; median={median_km:.3f} km"
    )

    haul_path = write_haul_csv(yearly, combined_mean_km, median_km, n)
    print(f"[OK] Haul table -> {haul_path}")

    for label, haul, filename in (
        ("mean", combined_mean_km, "figure11_mean.png"),
        ("median", median_km, "figure11_median.png"),
    ):
        china = _per_km_to_per_ton(*CHINA_TRANSPORT_PER_KM, haul)
        print(
            f"China transport ({label} {haul:.2f} km): "
            f"{china[0]:.2f}–{china[1]:.2f} USD/ton"
        )
        out = chart_path(filename)
        plot_figure11(out, data=component_data(haul))
        _export_high_res(out)


if __name__ == "__main__":
    main()
