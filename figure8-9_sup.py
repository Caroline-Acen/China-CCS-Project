"""
Figure 8-9 Supplementary: Cost Maps (Tank, avg cost type)
"""

from cost_maps import plot_cost_maps

plot_cost_maps(
    excel1_path="./data/DSA-Tank.xlsx",
    excel2_path="./data/EOR_Tank.xlsx",
    ne_admin0_dir="./data/natural_earth/admin0",
    ne_admin1_dir="./data/natural_earth/admin1",
    output_path="cost_maps_tank_avg_sup.png",
    cost_type="avg"
)
