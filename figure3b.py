"""
Figure 3b: Risk vs Storage Capacity Bubble Chart
Scatter plot showing Risk Factor vs Storage Capacity colored by storage type (DSA/EOR).
"""

import pandas as pd
import matplotlib.pyplot as plt


def plot_risk_vs_storage_bubble(
    csv_path,
    output_path="risk_capacity_tradeoff.png",
    dsa_color="green",
    eor_color="orange",
    figsize=(14, 9),
    s=55,
    alpha=0.75,
    dpi=300,
):
    """
    Create bubble chart: Risk Factor vs Storage Capacity.

    Parameters:
        csv_path: Path to CSV with storage_potential, final_risk_score, storage_type columns
        output_path: Output image path
        dsa_color: Color for DSA points
        eor_color: Color for EOR points
        figsize: Figure dimensions
        s: Marker size
        alpha: Marker transparency
        dpi: Output resolution
    """

    # Load and validate data
    df = pd.read_csv(csv_path, low_memory=False)
    required = {"storage_potential", "final_risk_score", "storage_type"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Clean data
    df = df[["storage_potential", "final_risk_score", "storage_type"]].copy()
    df["storage_potential"] = pd.to_numeric(df["storage_potential"], errors="coerce")
    df["final_risk_score"] = pd.to_numeric(df["final_risk_score"], errors="coerce")
    df["storage_type"] = df["storage_type"].astype(str).str.strip()
    df = df.dropna(subset=["storage_potential", "final_risk_score", "storage_type"])

    # Split by storage type
    dsa = df[df["storage_type"].str.upper() == "DSA"]
    eor = df[df["storage_type"].str.upper() == "EOR"]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    if not dsa.empty:
        ax.scatter(
            dsa["storage_potential"], 
            dsa["final_risk_score"],
            s=s, 
            c=dsa_color, 
            alpha=alpha,
            label="DSA"
        )

    if not eor.empty:
        ax.scatter(
            eor["storage_potential"], 
            eor["final_risk_score"],
            s=s, 
            c=eor_color, 
            alpha=alpha,
            label="EOR"
        )

    # Format axes
    ax.set_xlabel("Storage Capacity (Mt CO$_2$)", fontsize=18, fontweight="bold")
    ax.set_ylabel("Risk Factor", fontsize=18, fontweight="bold")
    ax.tick_params(axis="both", labelsize=12)
    ax.legend(loc="upper right", frameon=True, fontsize=13)
    ax.set_xlim(0, 350)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close(fig)
    print(f"[OK] Saved bubble chart to {output_path}")


# Run function
plot_risk_vs_storage_bubble(
    csv_path="./data/Risk_Assessment/final_seismic_risk_factor_base_case.csv",
    output_path="risk_storage_bubble.png"
)
