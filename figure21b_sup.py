"""Figure 21b Supplementary: Carbon Contributions Bar Charts (EOR)"""

from carbon_contributions import create_carbon_contribution_barcharts
from config import CARBON_EOR, chart_path

create_carbon_contribution_barcharts(
    file_path=str(CARBON_EOR),
    output_path=chart_path("carbon_contributions_eor.png"),
)
