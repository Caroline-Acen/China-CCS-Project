"""
Figure 3 (submitted): Cost difference maps — pipeline vs tank, average scenario.
"""

from config import DSA_PIPELINE, DSA_TANK, EOR_PIPELINE, EOR_TANK, NE_ADMIN0, NE_ADMIN1, chart_path
from cost_diff_maps import plot_cost_difference_maps

plot_cost_difference_maps(
    pipeline_excel1_path=str(DSA_PIPELINE),
    pipeline_excel2_path=str(EOR_PIPELINE),
    tank_excel1_path=str(DSA_TANK),
    tank_excel2_path=str(EOR_TANK),
    ne_admin0_dir=str(NE_ADMIN0),
    ne_admin1_dir=str(NE_ADMIN1),
    output_path=chart_path("Cost difference (EOR+DSA tank avg - EOR+DSA pipeline avg).png"),
    cost_type="avg",
    colorbar=2,
)
