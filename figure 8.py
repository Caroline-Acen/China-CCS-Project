import matplotlib.pyplot as plt
import numpy as np
import os

# Create directory if it doesn't exist
output_dir = "E:/CCS/Charts"
os.makedirs(output_dir, exist_ok=True)

# Data from the table
countries = ['China', 'USA', 'EU']
cost_components = ['Capture', 'Transportation', 'Storage']

# Cost ranges (USD/ton) - using the lower and upper bounds
cost_ranges = {
    'China': {
        'Capture': [32.02, 43.16],
        'Transportation': [0.11, 0.19],
        'Storage': [6.96, 8.35]
    },
    'USA': {
        'Capture': [25, 30],
        'Transportation': [25, 35],
        'Storage': [25, 35]
    },
    'EU': {
        'Capture': [40, 90],
        'Transportation': [2, 30],
        'Storage': [5, 35]
    }
}

# Set up the plot
fig, ax = plt.subplots(figsize=(14, 10))

# Colors for each cost component
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green for Capture, Transportation, Storage

# Create positions for each country
x_positions = np.arange(len(countries))

# Plot stacked bars for each country
for i, country in enumerate(countries):
    # Get cost ranges for this country
    capture_min, capture_max = cost_ranges[country]['Capture']
    transport_min, transport_max = cost_ranges[country]['Transportation']
    storage_min, storage_max = cost_ranges[country]['Storage']
    
    # Calculate heights for stacking
    capture_height = capture_max - capture_min
    transport_height = transport_max - transport_min
    storage_height = storage_max - storage_min
    
    # Bottom positions for stacking - each bar starts exactly where previous ends
    capture_bottom = capture_min
    transport_bottom = capture_max  # Starts exactly at top of capture bar
    storage_bottom = capture_max + transport_height  # Starts exactly at top of transport bar
    
    # Plot stacked bars with no gaps
    ax.bar(x_positions[i], capture_height, width=0.3, 
           bottom=capture_bottom, color=colors[0], alpha=0.7, 
           edgecolor='black', linewidth=1.5, label='Capture' if i == 0 else "")
    
    ax.bar(x_positions[i], transport_height, width=0.3, 
           bottom=transport_bottom, color=colors[1], alpha=0.7, 
           edgecolor='black', linewidth=1.5, label='Transportation' if i == 0 else "")
    
    ax.bar(x_positions[i], storage_height, width=0.3, 
           bottom=storage_bottom, color=colors[2], alpha=0.7, 
           edgecolor='black', linewidth=1.5, label='Storage' if i == 0 else "")

# Add vertical lines to divide each country section
for i in range(len(countries) - 1):
    # Calculate position between countries
    divider_x = i + 0.5
    ax.axvline(x=divider_x, color='gray', linestyle='--', alpha=0.7, linewidth=1)

# Customize the plot
ax.set_xlabel('Countries', fontsize=13, fontweight='bold')
ax.set_ylabel('Cost (USD/ton)', fontsize=13, fontweight='bold')
#ax.set_title('CCS Cost Stacked Comparison by Country', fontsize=14, fontweight='bold')

# Set x-ticks and labels
ax.set_xticks(x_positions)
ax.set_xticklabels(countries, fontsize=14)
#ax.set_xticklabels(countries, fontsize=16)

# Create custom legend
from matplotlib.patches import Patch

legend_elements = []
for component, color in zip(cost_components, colors):
    legend_elements.append(Patch(facecolor=color, alpha=0.7, edgecolor='black', label=component))

ax.legend(handles=legend_elements, loc='upper right', fontsize=16)

# Add grid for better readability
#ax.grid(True, alpha=0.3, axis='y')

# Adjust y-axis limit to accommodate all data
# Calculate maximum total cost for each country
max_total_costs = []
for country in countries:
    capture_max = cost_ranges[country]['Capture'][1]
    transport_max = cost_ranges[country]['Transportation'][1]
    storage_max = cost_ranges[country]['Storage'][1]
    total_max = capture_max + transport_max + storage_max
    max_total_costs.append(total_max)

max_value = max(max_total_costs)
min_value = min([cost_ranges[country]['Capture'][0] for country in countries])
ax.set_ylim(min_value * 0.9, max_value * 1.1)

# Adjust x-axis limits to show all countries clearly
ax.set_xlim(-0.5, len(countries) - 0.5)

# Add note about EU currency conversion
plt.figtext(0.02, 0.01, "Note: EU costs converted from EUR to USD: 1 USD = 0.87 EUR", 
           fontsize=14, style='italic')

plt.tight_layout()

# Save the figure
output_path = os.path.join(output_dir, "ccs_stacked_bar_comparison_by_country_4.png")
plt.savefig(output_path, dpi=600, bbox_inches='tight')
print(f"Figure saved to: {output_path}")

# Show the plot
plt.show()