"""Figure 21a Supplementary: Carbon Contributions Bar Charts (DSA)"""

from carbon_contributions import create_carbon_contribution_barcharts
from config import CARBON_DSA, chart_path

create_carbon_contribution_barcharts(
    file_path=str(CARBON_DSA),
    output_path=chart_path("carbon_contributions_dsa.png"),
)
