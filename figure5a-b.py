"""
Figure 5a-b: CO2 Flow Sankey Diagrams
Interactive Sankey diagrams showing CO2 flow allocation between provinces (DSA and EOR).
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.subplots as sp
import difflib


# Province bounding box coordinates (lat_min, lat_max, lon_min, lon_max)
PROVINCES = {
    'Beijing': (39, 41, 115, 117),
    'Tianjin': (38, 40, 116, 118),
    'Shanghai': (30, 32, 120, 122),
    'Chongqing': (28, 32, 105, 110),
    'Ningxia': (35, 39, 104, 107),
    'Hainan': (18, 21, 108, 111),
    'Hong Kong': (22, 23, 113, 115),
    'Macau': (22, 23, 112, 114),
    'Xinjiang': (35, 49, 75, 95),
    'Inner Mongolia': (37, 53, 97, 126),
    'Heilongjiang': (44, 54, 121, 135),
    'Jilin': (41, 46, 121, 131),
    'Liaoning': (38, 43, 118, 126),
    'Hebei': (36, 42, 113, 120),
    'Shandong': (34, 38, 114, 123),
    'Jiangsu': (31, 35, 116, 122),
    'Zhejiang': (27, 31, 118, 123),
    'Fujian': (23, 28, 115, 121),
    'Guangdong': (20, 25, 109, 117),
    'Sichuan': (26, 34, 97, 110),
    'Yunnan': (21, 29, 97, 106),
    'Tibet': (26, 36, 78, 99),
    'Qinghai': (31, 39, 89, 103),
    'Gansu': (32, 43, 92, 109),
    'Shaanxi': (31, 40, 105, 112),
    'Shanxi': (34, 41, 110, 115),
    'Henan': (31, 37, 110, 117),
    'Hubei': (29, 33, 108, 116),
    'Hunan': (24, 30, 108, 114),
    'Jiangxi': (24, 30, 113, 118),
    'Anhui': (29, 34, 114, 119),
    'Guangxi': (21, 26, 104, 112),
    'Guizhou': (24, 29, 103, 110)
}

# Province color mapping
PROVINCE_COLORS = {
    'Anhui': '#1f77b4', 'Beijing': '#ff7f0e', 'Fujian': '#2ca02c',
    'Gansu': '#d62728', 'Guangdong': '#9467bd', 'Guangxi': '#8c564b',
    'Guizhou': '#e377c2', 'Hainan': '#7f7f7f', 'Hebei': '#bcbd22',
    'Henan': '#17becf', 'Heilongjiang': '#aec7e8', 'Hubei': '#ffbb78',
    'Hunan': '#98df8a', 'Jilin': '#ff9896', 'Jiangsu': '#c5b0d5',
    'Jiangxi': '#c49c94', 'Liaoning': '#f7b6d2', 'Inner Mongolia': '#c7c7c7',
    'Ningxia': '#dbdb8d', 'Qinghai': '#9edae5', 'Shandong': '#393b79',
    'Shanxi': '#637939', 'Shaanxi': '#8c6d31', 'Shanghai': '#843c39',
    'Sichuan': '#7b4173', 'Tianjin': '#5254a3', 'Tibet': '#8ca252',
    'Xinjiang': '#bd9e39', 'Yunnan': '#ad494a', 'Zhejiang': '#a55194',
    'Chongqing': '#6b6ecf'
}


def normalize_province_name(name, provinces):
    """Normalize province names using fuzzy matching."""
    if pd.isna(name):
        return None
    clean = str(name).strip()
    if clean in provinces:
        return clean
    possible = difflib.get_close_matches(clean, provinces.keys(), n=1, cutoff=0.6)
    return possible[0] if possible else clean


def haversine_distance(prov_a, prov_b, provinces):
    """Calculate approximate distance between province centers."""
    prov_a = normalize_province_name(prov_a, provinces)
    prov_b = normalize_province_name(prov_b, provinces)
    if prov_a not in provinces or prov_b not in provinces:
        return 100
    lat_a = np.mean(provinces[prov_a][:2])
    lon_a = np.mean(provinces[prov_a][2:])
    lat_b = np.mean(provinces[prov_b][:2])
    lon_b = np.mean(provinces[prov_b][2:])
    return np.sqrt((lat_a - lat_b) ** 2 + (lon_a - lon_b) ** 2)


def allocate_co2_flow(emissions, storage, provinces):
    """Allocate CO2 flow from emitters to storage, prioritizing nearby locations."""
    flows = []
    storage_copy = storage.copy()
    
    for src, emit in emissions.items():
        remaining = emit
        sorted_sinks = sorted(storage_copy.keys(), key=lambda x: haversine_distance(src, x, provinces))
        for sink in sorted_sinks:
            if remaining <= 0:
                break
            if storage_copy[sink] <= 0:
                continue
            flow = min(remaining, storage_copy[sink])
            flows.append((src, sink, flow))
            storage_copy[sink] -= flow
            remaining -= flow
    return flows


def create_separate_sankey_plots(carbon_dsa_path, carbon_eor_path):
    """Create separate DSA and EOR Sankey plots with 6 subplots each."""

    # Read datasets
    dsa = pd.read_excel(carbon_dsa_path)
    eor = pd.read_excel(carbon_eor_path)

    # Clean province names
    dsa["Province name"] = dsa["Province name"].apply(lambda x: normalize_province_name(x, PROVINCES))
    eor["Province name"] = eor["Province name"].apply(lambda x: normalize_province_name(x, PROVINCES))

    years = ['2025', '2030', '2035', '2040', '2045', '2050']

    def create_sankey_plot(df, storage_type, years):
        """Create Sankey plot with 6 subplots for given storage type."""
        fig = sp.make_subplots(
            rows=2, cols=3,
            subplot_titles=years,
            vertical_spacing=0.15,
            horizontal_spacing=0.10,
            specs=[[{"type": "sankey"}, {"type": "sankey"}, {"type": "sankey"}],
                   [{"type": "sankey"}, {"type": "sankey"}, {"type": "sankey"}]]
        )

        # Collect all provinces for consistent coloring
        all_provinces = set()
        for year in years:
            emission_col = f'CO2_{year}'
            storage_col = f'P_{year}'
            mask = (df[storage_col] > 0) & (df[emission_col] > 0)
            plot_data = df[mask].copy()
            emissions = plot_data.set_index("Province name")[emission_col].to_dict()
            storage = plot_data.set_index("Province name")[storage_col].to_dict()
            all_provinces.update(emissions.keys())
            all_provinces.update(storage.keys())

        province_color_map = {
            province: PROVINCE_COLORS.get(province, '#808080') for province in all_provinces
        }

        # Generate Sankey diagram for each year
        for i, year in enumerate(years):
            row = i // 3 + 1
            col = i % 3 + 1

            emission_col = f'CO2_{year}'
            storage_col = f'P_{year}'
            mask = (df[storage_col] > 0) & (df[emission_col] > 0)
            plot_data = df[mask].copy()

            emissions = plot_data.set_index("Province name")[emission_col].to_dict()
            storage = plot_data.set_index("Province name")[storage_col].to_dict()
            flows = allocate_co2_flow(emissions, storage, PROVINCES)

            sources = [f[0] for f in flows]
            sinks = [f[1] for f in flows]
            values = [f[2] for f in flows]
            all_nodes = list(dict.fromkeys(sources + sinks))
            node_indices = {name: i for i, name in enumerate(all_nodes)}

            # Create link colors based on source province
            link_colors = []
            for source in sources:
                color = province_color_map[source]
                if color.startswith('#'):
                    color = color.lstrip('#')
                    rgb = tuple(int(color[i:i+2], 16) for i in (0, 2, 4))
                    link_colors.append(f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, 0.6)')
                else:
                    link_colors.append('rgba(100, 100, 100, 0.6)')

            sankey_trace = go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=30,
                    thickness=18,
                    line=dict(color=[province_color_map[node] for node in all_nodes], width=0.5),
                    label=all_nodes,
                    color=[province_color_map[node] for node in all_nodes]
                ),
                link=dict(
                    source=[node_indices[s] for s in sources],
                    target=[node_indices[t] for t in sinks],
                    value=values,
                    color=link_colors
                ),
                domain=dict(
                    x=[(col - 1) / 3 + 0.04, col / 3 - 0.04],
                    y=[(2 - row) / 2 + 0.05, (3 - row) / 2 - 0.05]
                )
            )

            fig.add_trace(sankey_trace, row=row, col=col)

        fig.update_layout(
            title_text=f"CO2 Flow - {storage_type}",
            font=dict(size=12),
            height=1000,
            width=1600,
            showlegend=False
        )

        for i, annotation in enumerate(fig.layout.annotations):
            annotation.text = years[i]
            annotation.font.size = 14

        return fig

    # Generate DSA Sankey plot
    print("Generating DSA Sankey plot...")
    dsa_fig = create_sankey_plot(dsa, "DSA", years)
    dsa_fig.write_html("DSA_Sankey_Plot.html", include_plotlyjs="cdn")
    print("[OK] DSA Sankey plot saved to: DSA_Sankey_Plot.html")

    # Generate EOR Sankey plot
    print("Generating EOR Sankey plot...")
    eor_fig = create_sankey_plot(eor, "EOR", years)
    eor_fig.write_html("EOR_Sankey_Plot.html", include_plotlyjs="cdn")
    print("[OK] EOR Sankey plot saved to: EOR_Sankey_Plot.html")

    return dsa_fig, eor_fig


# Run function
create_separate_sankey_plots(
    carbon_dsa_path="./data/Carbon_DSA.xlsx",
    carbon_eor_path="./data/Carbon_EOR.xlsx"
)
