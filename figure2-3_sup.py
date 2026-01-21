"""
Figure 2-3 Supplementary: Cost Maps (Pipeline, min cost type)
"""

from cost_maps import plot_cost_maps

plot_cost_maps(
    excel1_path="./data/DSA-Pipeline.xlsx",
    excel2_path="./data/EOR_Pipeline.xlsx",
    ne_admin0_dir="./data/natural_earth/admin0",
    ne_admin1_dir="./data/natural_earth/admin1",
    output_path="cost_maps_pipeline_min.png",
    cost_type="min"
)
