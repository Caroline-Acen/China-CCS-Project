"""
Figure 21b Supplementary: Carbon Contributions Bar Charts (EOR)
"""

from carbon_contributions import create_carbon_contribution_barcharts

create_carbon_contribution_barcharts(
    file_path="./data/Carbon_EOR.xlsx",
    output_path="carbon_contributions_eor.png"
)
