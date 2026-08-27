"""
Dynamic Injectivity-Distance Pareto Allocation (DIDPA)

Sink Suitability Score:
    SS_{i,j} = (R_{i,t} * S_{i,t}) / d_{i,j}^2

R and S are converted from Mt and Mt/year to tonnes and tonnes/year before scoring.
Matching uses a Dynamic Saturation and Marginal Disregard (SMD) loop.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree

from config import CARBON_DSA, CARBON_EOR, INDUSTRIES_CSV, NE_ADMIN1

EARTH_RADIUS_KM = 6371.0
MT_TO_T = 1_000_000.0

# Province centroid coordinates (lat_min, lat_max, lon_min, lon_max)
PROVINCE_CENTROIDS = {
    "Beijing": (39, 41, 115, 117),
    "Tianjin": (38, 40, 116, 118),
    "Shanghai": (30, 32, 120, 122),
    "Chongqing": (28, 32, 105, 110),
    "Ningxia": (35, 39, 104, 107),
    "Hainan": (18, 21, 108, 111),
    "Xinjiang": (35, 49, 75, 95),
    "Inner Mongolia": (37, 53, 97, 126),
    "Heilongjiang": (44, 54, 121, 135),
    "Jilin": (41, 46, 121, 131),
    "Liaoning": (38, 43, 118, 126),
    "Hebei": (36, 42, 113, 120),
    "Shandong": (34, 38, 114, 123),
    "Jiangsu": (31, 35, 116, 122),
    "Zhejiang": (27, 31, 118, 123),
    "Fujian": (23, 28, 115, 121),
    "Guangdong": (20, 25, 109, 117),
    "Sichuan": (26, 34, 97, 110),
    "Yunnan": (21, 29, 97, 106),
    "Tibet": (26, 36, 78, 99),
    "Qinghai": (31, 39, 89, 103),
    "Gansu": (32, 43, 92, 109),
    "Shaanxi": (31, 40, 105, 112),
    "Shanxi": (34, 41, 110, 115),
    "Henan": (31, 37, 110, 117),
    "Hubei": (29, 33, 108, 116),
    "Hunan": (24, 30, 108, 114),
    "Jiangxi": (24, 30, 113, 118),
    "Anhui": (29, 34, 114, 119),
    "Guangxi": (21, 26, 104, 112),
    "Guizhou": (24, 29, 103, 110),
}


@dataclass
class SinkState:
    r_remaining_t: float
    s_remaining_t: float
    province: str
    locked: bool = False


@dataclass
class AllocationResult:
    flows: list[tuple[str, str, float]]
    unallocated_t: float
    emitter_remaining: dict[int, float] = field(default_factory=dict)


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two WGS84 points."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return float(2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a)))


def sink_suitability_score(r_t_per_year: float, s_t: float, distance_km: float) -> float:
    if distance_km <= 0 or r_t_per_year <= 0 or s_t <= 0:
        return 0.0
    return (r_t_per_year * s_t) / (distance_km ** 2)


def _normalize_province(name: str | None) -> str | None:
    if name is None or (isinstance(name, float) and np.isnan(name)):
        return None
    clean = str(name).strip()
    aliases = {
        "Inner Mongolia Autonomous Region": "Inner Mongolia",
        "Inner Mongolia": "Inner Mongolia",
        "Nei Mongol": "Inner Mongolia",
        "Xinjiang Uygur Autonomous Region": "Xinjiang",
        "Xinjiang Uyghur Autonomous Region": "Xinjiang",
        "Xinjiang": "Xinjiang",
        "西藏自治区": "Tibet",
        "新疆维吾尔自治区": "Xinjiang",
        "内蒙古自治区": "Inner Mongolia",
        "广西壮族自治区": "Guangxi",
        "宁夏回族自治区": "Ningxia",
        "香港特别行政区": "Hong Kong",
        "澳门特别行政区": "Macau",
    }
    return aliases.get(clean, clean)


def province_from_bbox(x: float, y: float) -> str | None:
    """Fast province lookup using bounding boxes (used for sink aggregation)."""
    for province, (lat_min, lat_max, lon_min, lon_max) in PROVINCE_CENTROIDS.items():
        if lat_min <= y <= lat_max and lon_min <= x <= lon_max:
            return province
    return None


def assign_provinces_fast(points_df: pd.DataFrame) -> pd.Series:
    return points_df.apply(lambda row: province_from_bbox(row["X"], row["Y"]), axis=1)


def assign_provinces(points_df: pd.DataFrame, admin1_dir=NE_ADMIN1) -> pd.Series:
    """Spatially join lon/lat points to Natural Earth China provinces."""
    gdf = gpd.GeoDataFrame(
        points_df.copy(),
        geometry=gpd.points_from_xy(points_df["X"], points_df["Y"]),
        crs="EPSG:4326",
    )
    admin1 = gpd.read_file(f"{admin1_dir}/ne_50m_admin_1_states_provinces.shp")
    if "admin" in admin1.columns:
        china = admin1[admin1["admin"] == "China"].copy()
        name_col = "name"
    else:
        china = admin1[admin1["ADM0_A3"] == "CHN"].copy()
        name_col = "NAME" if "NAME" in china.columns else "name"

    joined = gpd.sjoin(gdf, china[[name_col, "geometry"]], how="left", predicate="within")
    return joined[name_col].map(_normalize_province)


def _storage_columns(year: int) -> tuple[str, str]:
    if year == 2025:
        return "Injection Rate", "Storage Potential"
    return f"R_{year}", f"S_{year} (t-1)"


def load_storage_sinks(storage_excel: str | Path, year: int, admin1_dir=NE_ADMIN1) -> pd.DataFrame:
    df = pd.read_excel(storage_excel)
    r_col, s_col = _storage_columns(year)
    for col in (r_col, s_col, "X", "Y"):
        if col not in df.columns:
            raise KeyError(f"Missing '{col}' in storage file for year {year}")

    sinks = df[["X", "Y", r_col, s_col]].copy()
    sinks = sinks[(sinks[r_col] > 0) & (sinks[s_col] > 0)].reset_index(drop=True)
    sinks["sink_id"] = np.arange(len(sinks))
    sinks["province"] = assign_provinces_fast(sinks)
    sinks["R_mt"] = sinks[r_col].astype(float)
    sinks["S_mt"] = sinks[s_col].astype(float)
    return sinks


def aggregate_sinks_by_province(sinks: pd.DataFrame) -> pd.DataFrame:
    """Sum injectivity and storage by province for fast Sankey-scale matching."""
    sinks = sinks.dropna(subset=["province"]).copy()
    grouped = sinks.groupby("province", as_index=False).agg(
        R_mt=("R_mt", "sum"),
        S_mt=("S_mt", "sum"),
        X=("X", "mean"),
        Y=("Y", "mean"),
    )
    grouped["sink_id"] = np.arange(len(grouped))
    return grouped


def load_province_emitters(carbon_path: str | Path, year: int) -> pd.DataFrame:
    """Province-centroid emitters using provincial CO2 totals (Mt/year)."""
    carbon = pd.read_excel(carbon_path)
    carbon["Province name"] = carbon["Province name"].str.strip()
    emission_col = f"CO2_{year}"
    if emission_col not in carbon.columns:
        raise KeyError(f"Missing {emission_col} in {carbon_path}")

    province_emissions = carbon.set_index("Province name")[emission_col].to_dict()
    rows = []
    for emitter_id, (province, (lat_min, lat_max, lon_min, lon_max)) in enumerate(PROVINCE_CENTROIDS.items()):
        e_mt = float(province_emissions.get(province, 0.0))
        if e_mt > 0:
            rows.append({
                "emitter_id": emitter_id,
                "province": province,
                "X": (lon_min + lon_max) / 2,
                "Y": (lat_min + lat_max) / 2,
                "E_mt": e_mt,
            })
    return pd.DataFrame(rows)


def downscale_provincial_emissions(
    carbon_path: str | Path,
    year: int,
    density_column: str = "SUM19",
    industries_csv: str | Path = INDUSTRIES_CSV,
    admin1_dir=NE_ADMIN1,
    active_only: bool = True,
) -> pd.DataFrame:
    """
    Assign emitter magnitudes E_{j,t} (Mt/year) to industry grid nodes by
    proportional downscaling of provincial CO2 totals using facility density.
    """
    carbon = pd.read_excel(carbon_path)
    carbon["Province name"] = carbon["Province name"].str.strip()
    emission_col = f"CO2_{year}"
    if emission_col not in carbon.columns:
        raise KeyError(f"Missing {emission_col} in {carbon_path}")

    province_emissions = carbon.set_index("Province name")[emission_col].to_dict()

    usecols = ["X", "Y", density_column]
    chunks = []
    for chunk in pd.read_csv(industries_csv, usecols=usecols, chunksize=250_000, low_memory=False):
        chunk = chunk[chunk[density_column] > 0] if active_only else chunk
        if not chunk.empty:
            chunks.append(chunk)
    if not chunks:
        raise ValueError("No industry grid cells found for emission downscaling.")

    industries = pd.concat(chunks, ignore_index=True)
    industries = (
        industries.groupby(["X", "Y"], as_index=False)[density_column]
        .sum()
        .rename(columns={density_column: "density"})
    )
    industries["province"] = assign_provinces(industries, admin1_dir=admin1_dir)
    industries = industries.dropna(subset=["province"])

    industries["prov_total_density"] = industries.groupby("province")["density"].transform("sum")
    industries["E_mt"] = industries.apply(
        lambda row: (
            province_emissions.get(row["province"], 0.0)
            * row["density"]
            / row["prov_total_density"]
            if row["prov_total_density"] > 0
            else 0.0
        ),
        axis=1,
    )
    industries = industries[industries["E_mt"] > 0].reset_index(drop=True)
    industries["emitter_id"] = np.arange(len(industries))
    return industries


def filter_emitters_near_storage(
    emitters: pd.DataFrame,
    sinks: pd.DataFrame,
    max_distance_km: float = 250.0,
) -> pd.DataFrame:
    """Keep only emitter grid cells within max_distance_km of at least one sink."""
    sink_coords = np.radians(sinks[["Y", "X"]].to_numpy())
    emitter_coords = np.radians(emitters[["Y", "X"]].to_numpy())
    tree = BallTree(sink_coords, metric="haversine")
    radius_rad = max_distance_km / EARTH_RADIUS_KM
    has_neighbor = tree.query_radius(emitter_coords, r=radius_rad, count_only=True) > 0
    filtered = emitters.loc[has_neighbor].copy().reset_index(drop=True)
    filtered["emitter_id"] = np.arange(len(filtered))
    return filtered


def run_didpa_smd(
    emitters: pd.DataFrame,
    sinks: pd.DataFrame,
    max_distance_km: float = 250.0,
) -> AllocationResult:
    """
    Run the SMD allocation loop.

    emitters: X, Y, E_mt, emitter_id, province
    sinks: X, Y, R_mt, S_mt, sink_id, province
    """
    emitter_ids = emitters["emitter_id"].to_numpy()
    emitter_remaining = {
        int(eid): float(emt) * MT_TO_T
        for eid, emt in zip(emitter_ids, emitters["E_mt"], strict=True)
    }

    sink_states: dict[int, SinkState] = {
        int(row.sink_id): SinkState(
            r_remaining_t=float(row.R_mt) * MT_TO_T,
            s_remaining_t=float(row.S_mt) * MT_TO_T,
            province=str(row.province) if pd.notna(row.province) else "Unknown",
        )
        for row in sinks.itertuples(index=False)
    }

    sink_coords = np.radians(sinks[["Y", "X"]].to_numpy())
    emitter_coords = np.radians(emitters[["Y", "X"]].to_numpy())
    tree = BallTree(sink_coords, metric="haversine")

    radius_rad = max_distance_km / EARTH_RADIUS_KM
    neighbor_lists = tree.query_radius(emitter_coords, r=radius_rad)

    # Precompute distances (km) for each emitter-sink neighbor pair
    neighbor_distances: list[list[tuple[int, float]]] = []
    for e_idx, sink_indices in enumerate(neighbor_lists):
        lat_e, lon_e = emitters.iloc[e_idx][["Y", "X"]]
        pairs = []
        for sink_idx in sink_indices:
            lat_s, lon_s = sinks.iloc[sink_idx][["Y", "X"]]
            pairs.append((int(sink_idx), haversine_km(lat_e, lon_e, lat_s, lon_s)))
        neighbor_distances.append(pairs)

    emitter_idx_by_id = {int(eid): idx for idx, eid in enumerate(emitter_ids)}
    pair_flows_t: dict[tuple[str, str], float] = {}

    while True:
        best_score = -1.0
        best_emitter = None
        best_sink = None

        for eid, e_rem in emitter_remaining.items():
            if e_rem <= 1e-6:
                continue
            e_idx = emitter_idx_by_id[eid]
            for sink_idx, dist_km in neighbor_distances[e_idx]:
                sink_id = int(sinks.iloc[sink_idx]["sink_id"])
                state = sink_states[sink_id]
                if state.locked:
                    continue
                score = sink_suitability_score(state.r_remaining_t, state.s_remaining_t, dist_km)
                if score > best_score:
                    best_score = score
                    best_emitter = eid
                    best_sink = sink_id

        if best_emitter is None or best_sink is None or best_score <= 0:
            break

        state = sink_states[best_sink]
        alloc_t = min(
            emitter_remaining[best_emitter],
            state.r_remaining_t,
            state.s_remaining_t,
        )
        if alloc_t <= 1e-6:
            state.locked = True
            continue

        src_prov = str(emitters.loc[emitters["emitter_id"] == best_emitter, "province"].iloc[0])
        pair = (src_prov, state.province)
        pair_flows_t[pair] = pair_flows_t.get(pair, 0.0) + alloc_t

        emitter_remaining[best_emitter] -= alloc_t
        state.r_remaining_t -= alloc_t
        state.s_remaining_t -= alloc_t

        if state.r_remaining_t <= 1e-6 or state.s_remaining_t <= 1e-6:
            state.locked = True

    flows = [(src, dst, amt / MT_TO_T) for (src, dst), amt in pair_flows_t.items()]
    unallocated_t = sum(max(v, 0.0) for v in emitter_remaining.values())
    return AllocationResult(
        flows=flows,
        unallocated_t=unallocated_t,
        emitter_remaining=emitter_remaining,
    )


def allocate_province_centroid_flows(
    storage_excel: str | Path,
    carbon_path: str | Path,
    year: int,
    max_distance_km: float = 250.0,
    admin1_dir=NE_ADMIN1,
) -> AllocationResult:
    """Fast DIDPA wrapper for Sankey diagrams (province emitters, province-aggregated sinks)."""
    emitters = load_province_emitters(carbon_path, year=year)
    sinks = aggregate_sinks_by_province(
        load_storage_sinks(storage_excel, year=year, admin1_dir=admin1_dir)
    )
    return run_didpa_smd(emitters, sinks, max_distance_km=max_distance_km)


def allocate_province_flows(
    storage_excel: str | Path,
    carbon_path: str | Path,
    year: int,
    density_column: str = "SUM19",
    max_distance_km: float = 250.0,
    admin1_dir=NE_ADMIN1,
) -> AllocationResult:
    """Convenience wrapper: downscale emissions, load sinks, run DIDPA."""
    emitters = downscale_provincial_emissions(
        carbon_path=carbon_path,
        year=year,
        density_column=density_column,
        admin1_dir=admin1_dir,
    )
    sinks = load_storage_sinks(storage_excel, year=year, admin1_dir=admin1_dir)
    emitters = filter_emitters_near_storage(emitters, sinks, max_distance_km=max_distance_km)
    return run_didpa_smd(emitters, sinks, max_distance_km=max_distance_km)
