"""
Figure 2 (submitted): Tank cost maps — DSA + EOR, average scenario.
"""

from config import DATA2_COST_UNIT, DSA_TANK, EOR_TANK, NE_ADMIN0, NE_ADMIN1, chart_path
from cost_maps import plot_cost_maps

plot_cost_maps(
    excel1_path=str(DSA_TANK),
    excel2_path=str(EOR_TANK),
    ne_admin0_dir=str(NE_ADMIN0),
    ne_admin1_dir=str(NE_ADMIN1),
    output_path=chart_path("EOR+DSA tank avg.png"),
    cost_type="avg",
    colorbar=2,
    apply_cny_conversion=False,
    cost_unit_label=DATA2_COST_UNIT,
)
