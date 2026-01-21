"""
Figure 4c: Storage Comparison (DSA vs EOR by Province)
Side-by-side scatter plots comparing storage capacity and injection rates by province.
"""

import pandas as pd
import matplotlib.pyplot as plt


def plot_storage_comparison(
    carbon_dsa_path, 
    carbon_eor_path, 
    output_path="storage_comparison.png"
):
    """
    Create two subplots comparing DSA and EOR storage potential and injection rates.

    Parameters:
        carbon_dsa_path: Path to Carbon DSA Excel file
        carbon_eor_path: Path to Carbon EOR Excel file
        output_path: Output image path
    """
    
    # Read data files
    carbon_dsa = pd.read_excel(carbon_dsa_path)
    carbon_eor = pd.read_excel(carbon_eor_path)
    
    # Get provinces with storage in either DSA or EOR
    all_provinces = set(carbon_dsa['Province name']).union(set(carbon_eor['Province name']))
    
    # Prepare plot data
    plot_data = []
    for province in all_provinces:
        # Get DSA data
        dsa_row = carbon_dsa[carbon_dsa['Province name'] == province]
        dsa_storage = dsa_row['Storage potential-DSA (Mt)'].values[0] if len(dsa_row) > 0 else 0
        dsa_injection = dsa_row['Injection rate capability-DSA (Mt/a) (Average)'].values[0] if len(dsa_row) > 0 else 0
        
        # Get EOR data
        eor_row = carbon_eor[carbon_eor['Province name'] == province]
        eor_storage = eor_row['Storage potential-EOR (Mt)'].values[0] if len(eor_row) > 0 else 0
        eor_injection = eor_row['Injection rate capability-EOR (Mt/a) (Average)'].values[0] if len(eor_row) > 0 else 0
        
        # Only include if province has storage
        if dsa_storage > 0 or eor_storage > 0:
            plot_data.append({
                'province': province,
                'dsa_storage': dsa_storage,
                'dsa_injection': dsa_injection,
                'eor_storage': eor_storage,
                'eor_injection': eor_injection
            })
    
    plot_df = pd.DataFrame(plot_data)
    
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
    
    # Create color dictionary with fallback
    color_dict = {
        province: province_colors.get(province.strip(), '#000000') 
        for province in plot_df['province']
    }
    
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), dpi=600)
    
    # Font sizes
    TITLE_FONT_SIZE = 21
    LABEL_FONT_SIZE = 18
    TICK_FONT_SIZE = 12
    LEGEND_FONT_SIZE = 15
    
    # DSA subplot
    dsa_data = plot_df[plot_df['dsa_storage'] > 0]
    for _, row in dsa_data.iterrows():
        ax1.scatter(
            row['dsa_storage'], row['dsa_injection'],
            c=[color_dict[row['province']]],
            marker='o', s=150, alpha=0.8,
            label=row['province']
        )
    
    ax1.set_xlabel('Storage Capacity (Mt CO$_2$)', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax1.set_ylabel('Injection Rate (Mt CO$_2$/year)', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax1.set_title('DSA', fontsize=TITLE_FONT_SIZE, fontweight='bold')
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(TICK_FONT_SIZE)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=LEGEND_FONT_SIZE)
    
    # EOR subplot
    eor_data = plot_df[plot_df['eor_storage'] > 0]
    for _, row in eor_data.iterrows():
        ax2.scatter(
            row['eor_storage'], row['eor_injection'],
            c=[color_dict[row['province']]],
            marker='s', s=150, alpha=0.8,
            label=row['province']
        )
    
    ax2.set_xlabel('Storage Capacity (Mt CO$_2$)', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax2.set_ylabel('Injection Rate (Mt CO$_2$/year)', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax2.set_title('EOR', fontsize=TITLE_FONT_SIZE, fontweight='bold')
    for label in ax2.get_xticklabels() + ax2.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(TICK_FONT_SIZE)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=LEGEND_FONT_SIZE)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()
    
    # Print summary
    print(f"[OK] Plot saved to {output_path}")
    print(f"Total provinces with DSA storage: {len(dsa_data)}")
    print(f"Total provinces with EOR storage: {len(eor_data)}")


# Run function
plot_storage_comparison(
    carbon_dsa_path="./data/Carbon_DSA.xlsx",
    carbon_eor_path="./data/Carbon_EOR.xlsx",
    output_path="storage_comparison.png"
)
