"""Supplementary Figure 6 and 7: Pipeline min-cost maps (DSA and EOR)."""

from config import chart_path
from supplementary_cost_maps import DATASETS, plot_dsa_cost_maps, plot_eor_cost_maps

plot_eor_cost_maps(
    eor_path=DATASETS["eor_pipeline"],
    dsa_path=DATASETS["dsa_pipeline"],
    output_path=chart_path("EOR pipeline min.png"),
    cost_type="min",
)

plot_dsa_cost_maps(
    dsa_path=DATASETS["dsa_pipeline"],
    output_path=chart_path("DSA pipeline min.png"),
    cost_type="min",
)
