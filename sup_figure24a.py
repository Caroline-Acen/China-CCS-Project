"""Supplementary Figure 24a: Carbon contribution bar charts (DSA)."""

from carbon_contributions import create_carbon_contribution_barcharts
from config import CARBON_DSA, chart_path

create_carbon_contribution_barcharts(
    file_path=str(CARBON_DSA),
    output_path=chart_path("sup_figure24a.png"),
)
