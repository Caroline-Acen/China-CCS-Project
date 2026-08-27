"""
Generate Supplementary Table 1: DSA full-chain cost statistics (pipeline transport).

Computes minimum, mean, maximum, and selected percentiles for each year and
cost scenario (min / avg / max) from the data grid files.
"""

from cost_maps import compute_cost_statistics
from config import DATA2_COST_UNIT, DSA_PIPELINE, data_path

if __name__ == "__main__":
    stats = compute_cost_statistics(
        DSA_PIPELINE,
        apply_cny_conversion=False,
        output_csv=data_path("supplementary_table1_dsa_pipeline.csv"),
        label=f"DSA pipeline ({DATA2_COST_UNIT})",
    )
    print(stats.to_string(index=False).replace("\u2082", "2"))
