"""
Figure 7 and Figure 8: CO2 flow Sankey diagrams — DSA and EOR.

Matches the approved left->right Plotly layout:
  - Years 2025-2050, 2x3 panels
  - Source provinces (left) -> Storage provinces (right)
  - Teal = Integrated Local Storage Hubs (source == sink)
  - Orange = Cross-Provincial Export Corridor (source != sink)
  - Province names on both sides
  - Storage Gt = cumulative stored since 2025 through that panel year
  - Ribbon hover = that year's flow (Gt); no mid-ribbon number boxes
"""

from __future__ import annotations

import difflib
from collections import defaultdict

import plotly.graph_objects as go
import plotly.subplots as sp

from config import (
    CARBON_DSA,
    CARBON_EOR,
    DSA_PIPELINE,
    EOR_PIPELINE,
    chart_path,
    manuscript_fontsize,
)
from didpa import allocate_province_centroid_flows
from plotly_export import write_plotly_html

# Okabe–Ito (colorblind-safe, journal-friendly)
COLOR_HUB = "rgba(0, 158, 115, 0.78)"       # #009E73 bluish green — hubs
COLOR_CORRIDOR = "rgba(230, 159, 0, 0.72)"  # #E69F00 orange — corridors
COLOR_NODE = "#9e9e9e"
COLOR_NODE_LINE = "#424242"
COLOR_HUB_LEGEND = "rgb(0, 158, 115)"
COLOR_CORRIDOR_LEGEND = "rgb(230, 159, 0)"

PROVINCE_NAMES = {
    "Anhui", "Beijing", "Fujian", "Gansu", "Guangdong", "Guangxi", "Guizhou",
    "Hainan", "Hebei", "Henan", "Heilongjiang", "Hubei", "Hunan", "Jilin",
    "Jiangsu", "Jiangxi", "Liaoning", "Inner Mongolia", "Ningxia", "Qinghai",
    "Shandong", "Shanxi", "Shaanxi", "Shanghai", "Sichuan", "Tianjin", "Tibet",
    "Xinjiang", "Yunnan", "Zhejiang", "Chongqing",
}

MT_TO_GT = 1_000.0
SANKEY_WIDTH_PX = 1700
SANKEY_HEIGHT_PX = 1900
SANKEY_WIDTH_IN = SANKEY_WIDTH_PX / 96.0
SANKEY_FS = round(manuscript_fontsize(SANKEY_WIDTH_IN), 1)
NODE_PAD = 36
NODE_THICKNESS = 14


def _fmt_stored_gt(gt: float) -> str:
    """Adaptive decimals so minute cumulative add-ups stay visible."""
    if gt <= 0:
        return "0 Gt"
    if gt >= 0.01:
        return f"{gt:.2f} Gt"
    if gt >= 0.0001:
        return f"{gt:.4f} Gt"
    return f"{gt:.6f}".rstrip("0").rstrip(".") + " Gt"


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


def _aggregate_flows(flows: list[tuple[str, str, float]]) -> list[tuple[str, str, float]]:
    totals: dict[tuple[str, str], float] = defaultdict(float)
    for src, sink, amt in flows:
        totals[(src, sink)] += amt
    return [(s, t, v) for (s, t), v in sorted(totals.items(), key=lambda kv: -kv[1])]


def _cumulative_stored_gt(
    year_flows: dict[int, list[tuple[str, str, float]]],
    years: list[int],
    upto_year: int,
) -> dict[str, float]:
    """Cumulative Gt stored in each province (as sink) from 2025 through upto_year."""
    stored_mt: dict[str, float] = defaultdict(float)
    for y in years:
        if y > upto_year:
            break
        for _, sink, amt_mt in year_flows[y]:
            stored_mt[sink] += amt_mt
    return {p: mt / MT_TO_GT for p, mt in stored_mt.items()}


def _build_year_sankey(
    flows_mt,
    cum_stored_gt: dict[str, float],
    domain: dict,
) -> go.Sankey:
    flows = _aggregate_flows(flows_mt)
    if not flows:
        return go.Sankey(
            node=dict(label=["No flow"], color=[COLOR_NODE]),
            link=dict(source=[], target=[], value=[]),
            domain=domain,
        )

    sources = [f[0] for f in flows]
    sinks = [f[1] for f in flows]
    values_gt = [f[2] / MT_TO_GT for f in flows]

    left_names = list(dict.fromkeys(sources))
    right_names = list(dict.fromkeys(sinks))
    n_left = len(left_names)

    left_labels = list(left_names)
    right_labels = [
        f"{p}  {_fmt_stored_gt(cum_stored_gt.get(p, 0.0))}" for p in right_names
    ]
    all_labels = left_labels + right_labels
    left_index = {p: i for i, p in enumerate(left_names)}
    right_index = {p: i + n_left for i, p in enumerate(right_names)}

    link_source = [left_index[s] for s, _, _ in flows]
    link_target = [right_index[t] for _, t, _ in flows]
    link_colors = [
        COLOR_HUB if s == t else COLOR_CORRIDOR for s, t, _ in flows
    ]

    return go.Sankey(
        arrangement="snap",
        orientation="h",
        valueformat=".3f",
        valuesuffix=" Gt",
        node=dict(
            pad=NODE_PAD,
            thickness=NODE_THICKNESS,
            line=dict(color=COLOR_NODE_LINE, width=0.5),
            label=all_labels,
            color=[COLOR_NODE] * len(all_labels),
            hoverinfo="skip",
        ),
        link=dict(
            source=link_source,
            target=link_target,
            value=values_gt,
            color=link_colors,
            hoverinfo="skip",
        ),
        domain=domain,
    )


def create_separate_sankey_plots(
    dsa_storage_path,
    eor_storage_path,
    carbon_dsa_path,
    carbon_eor_path,
    density_column="MC19",
):
    del density_column
    years = [2025, 2030, 2035, 2040, 2045, 2050]

    def create_sankey_plot(storage_excel, carbon_path, storage_type: str):
        print(f"  Computing DIDPA flows for {storage_type}...")
        year_flows: dict[int, list[tuple[str, str, float]]] = {}
        year_unallocated: dict[int, float] = {}
        for year in years:
            flows, unallocated = didpa_flows_for_year(storage_excel, carbon_path, year)
            year_flows[year] = flows
            year_unallocated[year] = unallocated
            print(f"    {year}: {len(flows)} pairs, unalloc={unallocated:.1f} Mt")

        fig = sp.make_subplots(
            rows=2,
            cols=3,
            subplot_titles=[str(y) for y in years],
            vertical_spacing=0.14,
            horizontal_spacing=0.08,
            specs=[[{"type": "sankey"}] * 3, [{"type": "sankey"}] * 3],
        )

        for i, year in enumerate(years):
            row = i // 3 + 1
            col = i % 3 + 1
            domain = dict(
                x=[(col - 1) / 3 + 0.04, col / 3 - 0.04],
                y=[0.52, 0.80] if row == 1 else [0.36, 0.48],
            )
            cum_gt = _cumulative_stored_gt(year_flows, years, year)
            trace = _build_year_sankey(year_flows[year], cum_gt, domain)
            fig.add_trace(trace, row=row, col=col)

        fig.update_layout(
            font=dict(size=SANKEY_FS),
            height=SANKEY_HEIGHT_PX,
            width=SANKEY_WIDTH_PX,
            showlegend=False,
            hovermode=False,
            margin=dict(t=150, b=390, l=40, r=40),
        )

        for i, annotation in enumerate(fig.layout.annotations):
            if i >= len(years):
                break
            row_idx = i // 3
            col_idx = i % 3
            year = years[i]
            unalloc_gt = year_unallocated[year] / MT_TO_GT
            annotation.text = f"{year}"
            if unalloc_gt > 0.0005:
                annotation.text += (
                    f"<br><span style='font-size:{SANKEY_FS * 0.85:.1f}px'>"
                    f"Unallocated: {unalloc_gt:.2f} Gt</span>"
                )
            annotation.font.size = SANKEY_FS
            annotation.font.color = "black"
            annotation.showarrow = False
            annotation.xref = "paper"
            annotation.yref = "paper"
            annotation.x = (col_idx + 0.5) / 3
            annotation.xanchor = "center"
            annotation.yanchor = "bottom"
            if row_idx == 0:
                annotation.y = 1.0
                annotation.yshift = 36
            else:
                annotation.y = 0.48
                annotation.yshift = 8

        fig.add_annotation(
            text=(
                f"<span style='color:{COLOR_HUB_LEGEND}'>■</span> Integrated Local Storage Hubs"
                " &nbsp;&nbsp; "
                f"<span style='color:{COLOR_CORRIDOR_LEGEND}'>■</span> Cross-Provincial Export Corridor"
            ),
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.0,
            yshift=-150,
            showarrow=False,
            font=dict(size=SANKEY_FS, color="#212121"),
            xanchor="center",
            align="center",
        )

        return fig

    print("Generating DSA Sankey plot...")
    dsa_fig = create_sankey_plot(str(DSA_PIPELINE), str(CARBON_DSA), "DSA")
    write_plotly_html(dsa_fig, chart_path("DSA_Sankey_Plot.html"))
    print(f"[OK] DSA Sankey plot saved to: {chart_path('DSA_Sankey_Plot.html')}")

    print("Generating EOR Sankey plot...")
    eor_fig = create_sankey_plot(str(EOR_PIPELINE), str(CARBON_EOR), "EOR")
    write_plotly_html(eor_fig, chart_path("EOR_Sankey_Plot.html"))
    print(f"[OK] EOR Sankey plot saved to: {chart_path('EOR_Sankey_Plot.html')}")

    return dsa_fig, eor_fig


if __name__ == "__main__":
    create_separate_sankey_plots(
        dsa_storage_path=str(DSA_PIPELINE),
        eor_storage_path=str(EOR_PIPELINE),
        carbon_dsa_path=str(CARBON_DSA),
        carbon_eor_path=str(CARBON_EOR),
        density_column="MC19",
    )
