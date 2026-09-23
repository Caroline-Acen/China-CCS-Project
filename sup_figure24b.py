"""Supplementary Figure 24b: Carbon contribution bar charts (EOR)."""

from carbon_contributions import create_carbon_contribution_barcharts
from config import CARBON_EOR, chart_path

create_carbon_contribution_barcharts(
    file_path=str(CARBON_EOR),
    output_path=chart_path("sup_figure24b.png"),
)
