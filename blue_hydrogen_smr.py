"""
Dual-Gate SMR optimization for blue hydrogen maps.

Gasfield-Gate (SMR at gas field g):
    D = d_{g,s} * w_CO2 + d_{g,c} * w_H2

Sink-Gate (SMR at storage sink s):
    D = d_{g,s} * w_Gas + d_{s,c} * w_H2

Distances are straight-line Euclidean on projected coordinates.
The 1.3 regional winding factor is deliberately omitted (Blue hydrogen redo.docx).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.spatial import cKDTree
from shapely.geometry import box

from config import (
    BLUE_H2_INDUSTRIES,
    BLUE_H2_RADIUS_M,
    CHINA_GAS_FIELDS,
    DSA_PIPELINE,
    EOR_PIPELINE,
    INDUSTRIES_CSV,
    W_CO2,
    W_GAS,
    W_H2,
)

TARGET_CRS = "ESRI:102012"


def _points_gdf(df: pd.DataFrame, x_col: str = "X", y_col: str = "Y", crs: str = "EPSG:4326") -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        df.copy(),
        geometry=gpd.points_from_xy(df[x_col], df[y_col]),
        crs=crs,
    )


def load_gas_facilities(gas_csv=CHINA_GAS_FIELDS) -> gpd.GeoDataFrame:
    gas_df = pd.read_csv(gas_csv, low_memory=False)
    facilities = (
        gas_df.groupby("Facility_N", as_index=False)
        .agg(X=("X", "mean"), Y=("Y", "mean"))
    )
    return _points_gdf(facilities)


def load_storage_sites(dsa_path=DSA_PIPELINE, eor_path=EOR_PIPELINE) -> gpd.GeoDataFrame:
    frames = []
    for path, storage_type in ((dsa_path, "DSA"), (eor_path, "EOR")):
        df = pd.read_excel(path)
        subset = df.loc[df["Storage Potential"] > 0, ["X", "Y", "Storage Potential"]].copy()
        subset["storage_type"] = storage_type
        frames.append(subset)
    storage = pd.concat(frames, ignore_index=True)
    return _points_gdf(storage)


def load_h2_consumers(
    industries_csv=INDUSTRIES_CSV,
    within_geometry=None,
    chunk_size: int = 200_000,
) -> gpd.GeoDataFrame:
    cols = ["X", "Y", *BLUE_H2_INDUSTRIES]
    chunks = []
    for chunk in pd.read_csv(industries_csv, usecols=cols, low_memory=False, chunksize=chunk_size):
        mask = (chunk[list(BLUE_H2_INDUSTRIES)] > 0).any(axis=1)
        filtered = chunk.loc[mask, cols]
        if filtered.empty:
            continue
        if within_geometry is not None:
            gdf = _points_gdf(filtered).to_crs(TARGET_CRS)
            filtered = filtered.loc[gdf.geometry.within(within_geometry).values]
        if not filtered.empty:
            chunks.append(filtered)

    if not chunks:
        return gpd.GeoDataFrame(columns=cols, geometry=[], crs="EPSG:4326")

    consumers = pd.concat(chunks, ignore_index=True)
    return _points_gdf(consumers)


def _build_tree(gdf: gpd.GeoDataFrame) -> tuple[cKDTree, np.ndarray]:
    coords = np.column_stack([gdf.geometry.x.to_numpy(), gdf.geometry.y.to_numpy()])
    return cKDTree(coords), coords


def _nearest_within(tree: cKDTree, coords: np.ndarray, point_xy, radius_m: float):
    dist, idx = tree.query(point_xy, k=1, distance_upper_bound=radius_m)
    if not np.isfinite(dist) or idx >= len(coords):
        return None, None
    return float(dist), int(idx)


def gasfield_gate_assignments(
    gas_gdf: gpd.GeoDataFrame,
    sink_gdf: gpd.GeoDataFrame,
    consumer_gdf: gpd.GeoDataFrame,
    w_co2: float = W_CO2,
    w_h2: float = W_H2,
    radius_m: float = BLUE_H2_RADIUS_M,
) -> pd.DataFrame:
    sink_tree, _ = _build_tree(sink_gdf)
    cons_tree, _ = _build_tree(consumer_gdf)
    sink_xy = np.column_stack([sink_gdf.geometry.x, sink_gdf.geometry.y])
    cons_xy = np.column_stack([consumer_gdf.geometry.x, consumer_gdf.geometry.y])
    rows = []

    for _, gas_row in gas_gdf.iterrows():
        xy = (gas_row.geometry.x, gas_row.geometry.y)
        d_gs, sink_idx = _nearest_within(sink_tree, sink_xy, xy, radius_m)
        d_gc, cons_idx = _nearest_within(cons_tree, cons_xy, xy, radius_m)
        if sink_idx is None or cons_idx is None:
            continue
        sink_row = sink_gdf.iloc[sink_idx]
        cons_row = consumer_gdf.iloc[cons_idx]
        rows.append(
            {
                "gate": "gasfield",
                "facility": gas_row["Facility_N"],
                "gas_x": gas_row.geometry.x,
                "gas_y": gas_row.geometry.y,
                "sink_x": sink_row.geometry.x,
                "sink_y": sink_row.geometry.y,
                "consumer_x": cons_row.geometry.x,
                "consumer_y": cons_row.geometry.y,
                "d_co2_km": d_gs / 1000.0,
                "d_h2_km": d_gc / 1000.0,
                "d_total_weighted": d_gs * w_co2 + d_gc * w_h2,
                "storage_type": sink_row["storage_type"],
                "storage_potential": float(sink_row["Storage Potential"]),
            }
        )

    return pd.DataFrame(rows)


def sink_gate_assignments(
    gas_gdf: gpd.GeoDataFrame,
    sink_gdf: gpd.GeoDataFrame,
    consumer_gdf: gpd.GeoDataFrame,
    w_gas: float = W_GAS,
    w_h2: float = W_H2,
    radius_m: float = BLUE_H2_RADIUS_M,
) -> pd.DataFrame:
    gas_tree, _ = _build_tree(gas_gdf)
    cons_tree, _ = _build_tree(consumer_gdf)
    gas_xy = np.column_stack([gas_gdf.geometry.x, gas_gdf.geometry.y])
    cons_xy = np.column_stack([consumer_gdf.geometry.x, consumer_gdf.geometry.y])
    rows = []

    for idx in range(len(sink_gdf)):
        sink_row = sink_gdf.iloc[idx]
        xy = (sink_row.geometry.x, sink_row.geometry.y)
        d_gs, gas_idx = _nearest_within(gas_tree, gas_xy, xy, radius_m)
        d_sc, cons_idx = _nearest_within(cons_tree, cons_xy, xy, radius_m)
        if gas_idx is None or cons_idx is None:
            continue
        gas_match = gas_gdf.iloc[gas_idx]
        cons_match = consumer_gdf.iloc[cons_idx]
        rows.append(
            {
                "gate": "sink",
                "facility": gas_match["Facility_N"],
                "gas_x": gas_match.geometry.x,
                "gas_y": gas_match.geometry.y,
                "sink_x": sink_row.geometry.x,
                "sink_y": sink_row.geometry.y,
                "consumer_x": cons_match.geometry.x,
                "consumer_y": cons_match.geometry.y,
                "d_gas_km": d_gs / 1000.0,
                "d_h2_km": d_sc / 1000.0,
                "d_total_weighted": d_gs * w_gas + d_sc * w_h2,
                "storage_type": sink_row["storage_type"],
                "storage_potential": sink_row["Storage Potential"],
            }
        )

    return pd.DataFrame(rows)


def prepare_blue_hydrogen_layers():
    from shapely.geometry import box

    gas_gdf = load_gas_facilities().to_crs(TARGET_CRS)
    sink_gdf = load_storage_sites().to_crs(TARGET_CRS)
    bounds = gpd.GeoSeries(
        pd.concat([gas_gdf.geometry, sink_gdf.geometry], ignore_index=True),
        crs=TARGET_CRS,
    ).total_bounds
    pad = BLUE_H2_RADIUS_M
    search_zone = box(bounds[0] - pad, bounds[1] - pad, bounds[2] + pad, bounds[3] + pad)
    consumer_gdf = load_h2_consumers(within_geometry=search_zone).to_crs(TARGET_CRS)
    return gas_gdf, sink_gdf, consumer_gdf
