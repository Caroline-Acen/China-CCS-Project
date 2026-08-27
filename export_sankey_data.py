"""
Export DIDPA province flow tables for the R Sankey figure.

Writes (under data/):
  sankey_flows_{dsa|eor}.csv
  sankey_cumulative_{dsa|eor}.csv
"""

from __future__ import annotations

import difflib
from collections import defaultdict
from pathlib import Path

import pandas as pd

from config import (
    CARBON_DSA,
    CARBON_EOR,
    DSA_PIPELINE,
    EOR_PIPELINE,
    data_path,
    ensure_data_dir,
)
from didpa import allocate_province_centroid_flows

PROVINCE_NAMES = {
    "Anhui", "Beijing", "Fujian", "Gansu", "Guangdong", "Guangxi", "Guizhou",
    "Hainan", "Hebei", "Henan", "Heilongjiang", "Hubei", "Hunan", "Jilin",
    "Jiangsu", "Jiangxi", "Liaoning", "Inner Mongolia", "Ningxia", "Qinghai",
    "Shandong", "Shanxi", "Shaanxi", "Shanghai", "Sichuan", "Tianjin", "Tibet",
    "Xinjiang", "Yunnan", "Zhejiang", "Chongqing",
}

MT_TO_GT = 1_000.0
YEARS = [2025, 2030, 2035, 2040, 2045, 2050]


def normalize_province_name(name):
    if name is None:
        return None
    clean = str(name).strip()
    aliases = {
        "Nei Mongol": "Inner Mongolia",
        "Xinjiang Uygur": "Xinjiang",
        "Xizang": "Tibet",
    }
    clean = aliases.get(clean, clean)
    if clean in PROVINCE_NAMES:
        return clean
    match = difflib.get_close_matches(clean, PROVINCE_NAMES, n=1, cutoff=0.6)
    return match[0] if match else clean


def didpa_flows_for_year(storage_excel, carbon_path, year):
    result = allocate_province_centroid_flows(
        storage_excel=storage_excel,
        carbon_path=carbon_path,
        year=year,
    )
    flows = []
    for src, sink, amount_mt in result.flows:
        src = normalize_province_name(src)
        sink = normalize_province_name(sink)
        if src and sink and amount_mt > 0:
            flows.append((src, sink, float(amount_mt)))
    unallocated_mt = result.unallocated_t / 1_000_000.0
    return flows, unallocated_mt


def _aggregate_flows(flows):
    totals = defaultdict(float)
    for src, sink, amt in flows:
        totals[(src, sink)] += amt
    return [(s, t, v) for (s, t), v in sorted(totals.items(), key=lambda kv: -kv[1])]


def export_scenario(storage_excel: str, carbon_path: str, tag: str) -> tuple[Path, Path]:
    flow_rows = []
    year_flows = {}
    for year in YEARS:
        flows, unalloc = didpa_flows_for_year(storage_excel, carbon_path, year)
        flows = _aggregate_flows(flows)
        year_flows[year] = flows
        print(f"  {tag} {year}: {len(flows)} pairs, unalloc={unalloc:.1f} Mt")
        for src, sink, amt_mt in flows:
            flow_rows.append(
                {
                    "scenario": tag.upper(),
                    "year": year,
                    "source": src,
                    "sink": sink,
                    "flow_mt": amt_mt,
                    "flow_gt": amt_mt / MT_TO_GT,
                    "flow_type": "hub" if src == sink else "corridor",
                    "unalloc_mt": unalloc,
                    "unalloc_gt": unalloc / MT_TO_GT,
                }
            )

    cum_rows = []
    stored = defaultdict(float)
    for year in YEARS:
        for _, sink, amt_mt in year_flows[year]:
            stored[sink] += amt_mt
        for province, amt_mt in sorted(stored.items()):
            cum_rows.append(
                {
                    "scenario": tag.upper(),
                    "year": year,
                    "province": province,
                    "cum_stored_mt": amt_mt,
                    "cum_stored_gt": amt_mt / MT_TO_GT,
                }
            )

    ensure_data_dir()
    flows_path = Path(data_path(f"sankey_flows_{tag}.csv"))
    cum_path = Path(data_path(f"sankey_cumulative_{tag}.csv"))
    pd.DataFrame(flow_rows).to_csv(flows_path, index=False)
    pd.DataFrame(cum_rows).to_csv(cum_path, index=False)
    print(f"[OK] Wrote {flows_path}")
    print(f"[OK] Wrote {cum_path}")
    return flows_path, cum_path


def main():
    print("Exporting DSA Sankey tables...")
    export_scenario(str(DSA_PIPELINE), str(CARBON_DSA), "dsa")
    print("Exporting EOR Sankey tables...")
    export_scenario(str(EOR_PIPELINE), str(CARBON_EOR), "eor")


if __name__ == "__main__":
    main()
