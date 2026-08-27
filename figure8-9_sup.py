"""Figure 8-9 Supplementary: Tank avg-cost maps (data)."""

from config import chart_path
from supplementary_cost_maps import DATASETS, plot_dsa_cost_maps, plot_eor_cost_maps

plot_eor_cost_maps(
    eor_path=DATASETS["eor_tank"],
    dsa_path=DATASETS["dsa_tank"],
    output_path=chart_path("EOR tank avg.png"),
    cost_type="avg",
)

plot_dsa_cost_maps(
    dsa_path=DATASETS["dsa_tank"],
    output_path=chart_path("DSA tank avg.png"),
    cost_type="avg",
)
