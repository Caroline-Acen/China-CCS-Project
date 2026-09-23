"""Supplementary Figure 8 and 9: Tank min-cost maps (DSA and EOR)."""

from config import chart_path
from supplementary_cost_maps import DATASETS, plot_dsa_cost_maps, plot_eor_cost_maps

plot_eor_cost_maps(
    eor_path=DATASETS["eor_tank"],
    dsa_path=DATASETS["dsa_tank"],
    output_path=chart_path("EOR tank min.png"),
    cost_type="min",
)

plot_dsa_cost_maps(
    dsa_path=DATASETS["dsa_tank"],
    output_path=chart_path("DSA tank min.png"),
    cost_type="min",
)
