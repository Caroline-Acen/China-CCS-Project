"""
Figure 4d: Injection vs Unabated Emissions
Bar chart showing emissions and injection rates over time.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import ScalarFormatter


def plot_injection_unabated(
    data_file_path,
    output_path="injection_unabated.png"
):
    """
    Create side-by-side bar chart with two y-axes.

    Parameters:
        data_file_path: Path to Excel file with Year, Injection, Unabated Emissions columns
        output_path: Output image path

    Y-axes:
        - Left (red): Unabated Emissions
        - Right (blue): Injection Rate
        - Black line: Total Emissions = Unabated + Injection
    """
    
    # Read data
    df = pd.read_excel(data_file_path)
    
    # Validate columns
    required_columns = ['Year', 'Injection', 'Unabated Emissions']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    
    df = df.sort_values('Year')
    
    # Setup figure
    fig, ax1 = plt.subplots(figsize=(12, 8))
    years = df['Year'].values
    x_pos = np.arange(len(years))
    bar_width = 0.35
    
    # Secondary axis
    ax2 = ax1.twinx()
    
    # Unabated Emissions bars (left axis)
    ax1.bar(
        x_pos - bar_width/2, 
        df['Unabated Emissions'], 
        bar_width,
        label='Unabated Emissions', 
        color='red', 
        alpha=0.8
    )
    
    # Scale factor for visual balance
    scale_factor = 6
    scaled_injection = df['Injection'] * scale_factor
    
    # Injection Rate bars (right axis)
    ax2.bar(
        x_pos + bar_width/2, 
        scaled_injection, 
        bar_width,
        label='Injection Rate', 
        color='blue', 
        alpha=0.8
    )
    
    # Total Emissions line (smoothed)
    total_emissions = df['Unabated Emissions'] + df['Injection']
    x_new = np.linspace(x_pos.min(), x_pos.max(), 300)
    spl = make_interp_spline(x_pos, total_emissions, k=3)
    y_smooth = spl(x_new)
    ax1.plot(x_new, y_smooth, label='Total Emissions', color='black', linewidth=2)
    
    # Format axes
    ax1.set_xlabel('Year', fontsize=12)
    ax1.set_ylabel('Unabated Emissions / Total Emissions (Mt CO$_2$ / year)', color='red', fontsize=12)
    ax2.set_ylabel('Injection Rate (Mt CO$_2$ / year)', color='blue', fontsize=12)
    
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([str(int(year)) for year in years])
    
    # Scientific notation for y-axes
    for ax, color in [(ax1, 'red'), (ax2, 'blue')]:
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        ax.yaxis.set_major_formatter(formatter)
        ax.tick_params(axis='y', labelcolor=color)
    
    ax2.set_ylim(0, max(scaled_injection) * 2.0)
    
    # Combined legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"[OK] Plot saved to {output_path}")


# Run function
plot_injection_unabated(
    data_file_path="./data/Emission_Storage.xlsx",
    output_path="injection_unabated.png"
)
