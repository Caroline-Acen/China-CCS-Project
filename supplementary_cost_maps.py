"""Shared helpers for supplementary cost-map figure scripts."""

from config import (
    DATA2_COST_UNIT,
    DSA_PIPELINE,
    DSA_TANK,
    EOR_PIPELINE,
    EOR_TANK,
    NE_ADMIN0,
    NE_ADMIN1,
    chart_path,
)
from cost_maps import plot_cost_maps

COMMON_KWARGS = dict(
    ne_admin0_dir=str(NE_ADMIN0),
    ne_admin1_dir=str(NE_ADMIN1),
    colorbar=2,
    apply_cny_conversion=False,
    cost_unit_label=DATA2_COST_UNIT,
)

DATASETS = {
    "dsa_pipeline": str(DSA_PIPELINE),
    "eor_pipeline": str(EOR_PIPELINE),
    "dsa_tank": str(DSA_TANK),
    "eor_tank": str(EOR_TANK),
}


def plot_eor_cost_maps(*, eor_path: str, dsa_path: str, output_path: str, cost_type: str):
    """EOR-only cost maps that keep DSA area outlines where EOR fill is empty/NaN."""
    plot_cost_maps(
        excel1_path=eor_path,
        excel2_path=eor_path,
        output_path=output_path,
        cost_type=cost_type,
        outline_coverage_paths=(dsa_path, eor_path),
        **COMMON_KWARGS,
    )


def plot_dsa_cost_maps(*, dsa_path: str, output_path: str, cost_type: str):
    plot_cost_maps(
        excel1_path=dsa_path,
        excel2_path=dsa_path,
        output_path=output_path,
        cost_type=cost_type,
        **COMMON_KWARGS,
    )
