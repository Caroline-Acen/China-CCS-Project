"""Shared province color/marker styling for multi-province figures."""

import pandas as pd
import matplotlib.pyplot as plt

from config import CARBON_DSA, CARBON_EOR

MARKERS = ["o", "s", "^", "D", "v", "P", "*", "X", "h", "8", "<", ">", "p", "H", "d"]
STRONG_COLORS = plt.cm.tab20.colors


def build_province_style(provinces):
    sorted_provinces = sorted({str(p).strip() for p in provinces if str(p).strip()})
    return {
        p: {
            "color": STRONG_COLORS[i % len(STRONG_COLORS)],
            "marker": MARKERS[i % len(MARKERS)],
        }
        for i, p in enumerate(sorted_provinces)
    }


def province_style_from_carbon(
    carbon_dsa_path=CARBON_DSA,
    carbon_eor_path=CARBON_EOR,
):
    """Province styling aligned with figure4c storage comparison charts."""
    carbon_dsa = pd.read_excel(carbon_dsa_path)
    carbon_eor = pd.read_excel(carbon_eor_path)
    all_provinces = sorted(set(carbon_dsa["Province name"]).union(set(carbon_eor["Province name"])))
    style_provinces = []
    for province in all_provinces:
        dsa_row = carbon_dsa[carbon_dsa["Province name"] == province]
        eor_row = carbon_eor[carbon_eor["Province name"] == province]
        dsa_storage = dsa_row["Storage potential-DSA (Mt)"].values[0] if len(dsa_row) else 0
        eor_storage = eor_row["Storage potential-EOR (Mt)"].values[0] if len(eor_row) else 0
        if dsa_storage > 0 or eor_storage > 0:
            style_provinces.append(province.strip())
    return build_province_style(style_provinces)
