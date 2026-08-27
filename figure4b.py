"""
Figure 5b (submitted): Storage vs Injection Rate scatter — DSA and EOR.
Reviewer fix: fitted regression line to support stated trend.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import DSA_PIPELINE, EOR_PIPELINE, chart_path, manuscript_font_bundle


def plot_storage_vs_injection(
    dsa_excel_path,
    eor_excel_path,
    output_path="storage_vs_injection.png",
    dpi=300,
):
    df_dsa = pd.read_excel(dsa_excel_path)
    df_eor = pd.read_excel(eor_excel_path)

    storage_col = "Storage Potential"
    injection_col = "Injection Rate"

    fig, ax = plt.subplots(figsize=(10, 7))
    mf = manuscript_font_bundle(fig.get_figwidth())

    for df, color, label in [(df_eor, "orange", "EOR"), (df_dsa, "green", "DSA")]:
        mask = (df[storage_col] > 0) & (df[injection_col] > 0)
        subset = df.loc[mask]
        ax.scatter(
            subset[storage_col],
            subset[injection_col],
            c=color,
            s=30,
            alpha=0.55,
            label=label,
            edgecolors="none",
        )

        if len(subset) >= 2:
            x = subset[storage_col].to_numpy()
            y = subset[injection_col].to_numpy()
            slope, intercept = np.polyfit(x, y, 1)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax.plot(
                x_line,
                slope * x_line + intercept,
                color=color,
                linewidth=2,
                linestyle="--",
                alpha=0.9,
                label=f"{label} linear fit",
            )

    ax.set_xlabel("Total Storage Capacity (Mt CO$_2$)", fontsize=mf["label"])
    ax.set_ylabel("Injection Rate Capacity (Mt CO$_2$ / year)", fontsize=mf["label"])
    ax.set_xlim(0, 350)
    ax.set_ylim(0, 325)
    ax.legend(loc="upper right", fontsize=mf["legend"])
    ax.tick_params(axis="both", labelsize=mf["tick"])

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved scatter plot to {output_path}")


if __name__ == "__main__":
    plot_storage_vs_injection(
        dsa_excel_path=str(DSA_PIPELINE),
        eor_excel_path=str(EOR_PIPELINE),
        output_path=chart_path("storage_vs_injection.png"),
    )
