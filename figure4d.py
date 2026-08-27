"""
Figure 5d (submitted): Injection vs unabated emissions pathway.
Reviewer fix: single shared y-axis scale and larger labels.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import FuncFormatter

from config import EMISSION_STORAGE, chart_path, manuscript_font_bundle


def plot_injection_unabated(data_file_path, output_path="injection_unabated.png"):
    df = pd.read_excel(data_file_path)
    required = ["Year", "Injection", "Unabated Emissions"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df.sort_values("Year")
    total_emissions = df["Unabated Emissions"] + df["Injection"]

    fig, ax = plt.subplots(figsize=(12, 8))
    mf = manuscript_font_bundle(fig.get_figwidth())
    ax2 = ax.twinx()
    years = df["Year"].values
    x_pos = np.arange(len(years))
    width = 0.38

    unabated_bar = ax.bar(
        x_pos - width / 2,
        df["Unabated Emissions"],
        width,
        label="Unabated Emissions",
        color="#d62728",
        alpha=0.85,
    )
    injection_bar = ax2.bar(
        x_pos + width / 2,
        df["Injection"],
        width,
        label="Injection Rate",
        color="#0057b7",
        alpha=0.85,
    )

    spl = make_interp_spline(x_pos, total_emissions, k=min(3, len(x_pos) - 1))
    x_smooth = np.linspace(x_pos.min(), x_pos.max(), 200)
    total_line = ax.plot(
        x_smooth,
        spl(x_smooth),
        color="black",
        linewidth=2.5,
        label="Total Emissions",
    )

    ax.set_xlabel("Year", fontsize=mf["label"])
    ax.set_ylabel("Unabated emissions / total emissions (Mt CO$_2$/year)", fontsize=mf["label"], color="#d62728")
    ax2.set_ylabel("Captured / injected (Mt CO$_2$/year)", fontsize=mf["label"], color="#0057b7")
    ax.set_xticks(x_pos)
    ax.set_xticklabels([str(int(y)) for y in years], fontsize=mf["tick"])
    ax.tick_params(axis="y", labelsize=mf["tick"], colors="#d62728")
    ax2.tick_params(axis="y", labelsize=mf["tick"], colors="#0057b7")
    ax.set_ylim(0, max(total_emissions.max(), df["Unabated Emissions"].max()) * 1.15)
    ax2.set_ylim(0, df["Injection"].max() * 1.25)

    # Force plain numeric formatting (no 10^3 / 10^4 scientific scaling).
    comma_fmt = FuncFormatter(lambda x, _pos: f"{int(round(x)):,}")
    ax.yaxis.set_major_formatter(comma_fmt)
    ax2.yaxis.set_major_formatter(comma_fmt)

    # Exact legend entries: Total Emissions, Unabated Emissions, Injection Rate.
    ax.legend(
        [total_line[0], unabated_bar, injection_bar],
        ["Total Emissions", "Unabated Emissions", "Injection Rate"],
        loc="upper right",
        fontsize=mf["legend"],
        frameon=True,
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"[OK] Plot saved to {output_path}")


if __name__ == "__main__":
    plot_injection_unabated(
        data_file_path=str(EMISSION_STORAGE),
        output_path=chart_path("injection_unabated.png"),
    )
