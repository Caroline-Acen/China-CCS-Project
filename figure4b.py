"""
Figure 4b: Storage vs Injection Rate Scatter Plot
Scatter plot comparing Total Storage Capacity vs Injection Rate Capacity for DSA and EOR.
"""

import pandas as pd
import matplotlib.pyplot as plt


def plot_storage_vs_injection(
    dsa_excel_path,
    eor_excel_path,
    output_path="storage_vs_injection.png",
    dpi=300,
):
    """
    Create scatter plot: Total Storage Capacity vs Injection Rate Capacity.

    Parameters:
        dsa_excel_path: Path to DSA Excel file
        eor_excel_path: Path to EOR Excel file
        output_path: Output image path
        dpi: Output resolution
    """
    
    # Load data
    df_dsa = pd.read_excel(dsa_excel_path)
    df_eor = pd.read_excel(eor_excel_path)
    
    # Column names
    storage_col = "Storage Potential"
    injection_col = "Injection Rate"
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot EOR first (orange) so DSA (green) appears on top
    ax.scatter(
        df_eor[storage_col],
        df_eor[injection_col],
        c="orange",
        s=30,
        alpha=0.7,
        label="EOR",
        edgecolors="none"
    )
    
    # Plot DSA (green)
    ax.scatter(
        df_dsa[storage_col],
        df_dsa[injection_col],
        c="green",
        s=30,
        alpha=0.7,
        label="DSA",
        edgecolors="none"
    )
    
    # Format axes
    ax.set_xlabel("Total Storage Capacity (Mt CO$_2$)", fontsize=12)
    ax.set_ylabel("Injection Rate Capacity (Mt CO$_2$ / year)", fontsize=12)
    ax.set_xlim(0, 350)
    ax.set_ylim(0, 325)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(False)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved scatter plot to {output_path}")


# Run function
plot_storage_vs_injection(
    dsa_excel_path="./data/DSA-Pipeline.xlsx",
    eor_excel_path="./data/EOR_Pipeline.xlsx",
    output_path="storage_vs_injection.png"
)
