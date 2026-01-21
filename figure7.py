"""
Figure 7: Carbon Bubble Charts
Bubble charts showing CO2 emissions vs storage potential by province and year for DSA and EOR.
"""

import pandas as pd
import matplotlib.pyplot as plt


def generate_bubble_charts(dsa_file, eor_file, output_dir='.'):
    """
    Generate bubble charts for DSA and EOR data.
    
    Parameters:
        dsa_file: Path to Carbon_DSA.xlsx file
        eor_file: Path to Carbon_EOR.xlsx file
        output_dir: Directory to save output PNG files
    """
    
    # Read Excel files
    dsa_df = pd.read_excel(dsa_file)
    eor_df = pd.read_excel(eor_file)
    
    # Clean province names
    dsa_df['Province name'] = dsa_df['Province name'].str.strip()
    eor_df['Province name'] = eor_df['Province name'].str.strip()
    
    years = ['2025', '2030', '2035', '2040', '2045', '2050']
    
    # Province color mapping
    province_colors = {
        'Anhui': '#1f77b4', 'Beijing': '#ff7f0e', 'Fujian': '#2ca02c',
        'Gansu': '#d62728', 'Guangdong': '#9467bd', 'Guangxi': '#8c564b',
        'Guizhou': '#e377c2', 'Hainan': '#7f7f7f', 'Hebei': '#bcbd22',
        'Henan': '#17becf', 'Heilongjiang': '#aec7e8', 'Hubei': '#ffbb78',
        'Hunan': '#98df8a', 'Jilin': '#ff9896', 'Jiangsu': '#c5b0d5',
        'Jiangxi': '#c49c94', 'Liaoning': '#f7b6d2', 'Inner Mongolia': '#c7c7c7',
        'Ningxia': '#dbdb8d', 'Qinghai': '#9edae5', 'Shandong': '#393b79',
        'Shanxi': '#637939', 'Shaanxi': '#8c6d31', 'Shanghai': '#843c39',
        'Sichuan': '#7b4173', 'Tianjin': '#5254a3', 'Tibet': '#8ca252',
        'Xinjiang': '#bd9e39', 'Yunnan': '#ad494a', 'Zhejiang': '#a55194',
        'Chongqing': '#6b6ecf'
    }
    
    def create_bubble_plot(df, storage_type, years, province_colors):
        """Create bubble plot for given dataframe and storage type."""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        # Calculate bubble sizes based on injection rate
        injection_col = f'Injection rate capability-{storage_type} (Mt/a) (Average)'
        injection_rates = df[injection_col].values
        
        # Normalize for bubble sizes
        rate_range = injection_rates.max() - injection_rates.min()
        if rate_range > 0:
            bubble_sizes = (injection_rates - injection_rates.min()) / rate_range * 500 + 50
        else:
            bubble_sizes = [100] * len(injection_rates)
        
        for i, year in enumerate(years):
            ax = axes[i]
            emission_col = f'CO2_{year}'
            storage_col = f'P_{year}'

            mask = (df[storage_col] > 0) & (df[emission_col] > 0)
            plot_data = df[mask].copy()

            for province in plot_data['Province name']:
                province_data = plot_data[plot_data['Province name'] == province]
                if len(province_data) > 0:
                    emission = province_data[emission_col].values[0]
                    storage = province_data[storage_col].values[0]
                    idx = province_data.index[0]
                    
                    ax.scatter(
                        emission, storage, 
                        s=bubble_sizes[idx], 
                        c=[province_colors.get(province, '#000000')], 
                        alpha=0.7, 
                        label=province if i == 0 else ""
                    )

            ax.set_xlabel('CO2 Emissions (Mt)', fontsize=12)
            ax.set_ylabel('Storage Potential (Mt)', fontsize=12)
            ax.set_title(year, fontsize=14, fontweight='bold')

            # Use log scale if data spans multiple orders of magnitude
            if len(plot_data) > 0:
                if plot_data[emission_col].max() / max(plot_data[emission_col].min(), 1e-6) > 100:
                    ax.set_xscale('log')
                if plot_data[storage_col].max() / max(plot_data[storage_col].min(), 1e-6) > 100:
                    ax.set_yscale('log')

            # Add diagonal reference line
            ax.autoscale(enable=False)
            x_min, x_max = ax.get_xlim()
            y_min, y_max = ax.get_ylim()
            ax.plot([x_min, x_max], [y_min, y_max], color='red', alpha=0.2, linewidth=2, zorder=1)
        
        # Create legend
        handles, labels = axes[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        
        uniform_handles = [
            plt.Line2D([0], [0], marker='o', color=province_colors.get(label, '#000000'), 
                      markersize=10, linestyle='', alpha=0.7)
            for label in by_label.keys()
        ]
        
        n_cols = (len(by_label) + 3) // 4
        fig.legend(
            uniform_handles, by_label.keys(), 
            loc='lower center', 
            bbox_to_anchor=(0.5, 0.02), 
            ncol=n_cols, 
            fontsize=10,
            frameon=True,
            fancybox=True,
            shadow=True
        )
        
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.22, hspace=0.3)
        return fig
    
    # Generate DSA plot
    print("Generating DSA bubble chart...")
    dsa_fig = create_bubble_plot(dsa_df, 'DSA', years, province_colors)
    dsa_output_path = f'{output_dir}/Carbon_DSA_bubble_chart.png'
    dsa_fig.savefig(dsa_output_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(dsa_fig)
    print(f"[OK] DSA chart saved to: {dsa_output_path}")
    
    # Generate EOR plot
    print("Generating EOR bubble chart...")
    eor_fig = create_bubble_plot(eor_df, 'EOR', years, province_colors)
    eor_output_path = f'{output_dir}/Carbon_EOR_bubble_chart.png'
    eor_fig.savefig(eor_output_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(eor_fig)
    print(f"[OK] EOR chart saved to: {eor_output_path}")
    
    print("All charts generated successfully!")
    

# Run function
generate_bubble_charts(
    dsa_file='./data/Carbon_DSA.xlsx', 
    eor_file='./data/Carbon_EOR.xlsx'
)
