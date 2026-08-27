"""
Figure 1 (submitted): Pipeline cost maps — DSA + EOR, average scenario.
Also writes Supplementary Table 1 statistics for DSA pipeline.
"""

from config import DSA_PIPELINE, DATA2_COST_UNIT, EOR_PIPELINE, NE_ADMIN0, NE_ADMIN1, chart_path, data_path
from cost_maps import plot_cost_maps

plot_cost_maps(
    excel1_path=str(DSA_PIPELINE),
    excel2_path=str(EOR_PIPELINE),
    ne_admin0_dir=str(NE_ADMIN0),
    ne_admin1_dir=str(NE_ADMIN1),
    output_path=chart_path("EOR+DSA pipeline avg.png"),
    cost_type="avg",
    colorbar=2,
    apply_cny_conversion=False,
    cost_unit_label=DATA2_COST_UNIT,
    stats_csv_path=data_path("supplementary_table1_dsa_pipeline.csv"),
    stats_label="DSA pipeline",
)
