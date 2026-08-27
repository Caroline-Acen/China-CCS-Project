"""
Seismic risk factor R (Equation 5, EE-ART manuscript) with updated directional SFP.

    R = w1*SPGA + w2*SFP + w3*SSA     (w1,w2,w3) = (0.5, 0.3, 0.2)

Updated fault proximity (Risk Assessment Redo.docx):
    SFP_raw = FAIS * exp(-d/λ) * (1 + cos(α)) / 2
    SFP     = SFP_raw / FAIS_max * 10,   FAIS_max = 10

    FAIS = W_age * 1.5 + W_fea           (Table 3 / Eq. 7)
    α    = acute angle (0–90°) between site→fault vector and fault strike
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points
from shapely.strtree import STRtree
from sklearn.neighbors import BallTree

from config import (
    CAFD_FAULTS_CSV,
    DSA_PIPELINE,
    EOR_PIPELINE,
    PGA_CSV_BY_LEVEL,
    VS30_CSV,
)

CHINA_ALBERS = "EPSG:4479"
W1, W2, W3 = 0.5, 0.3, 0.2
FAIS_MAX = 10.0
BOORE_B = -1.186
BOORE_VREF = 760.0


@dataclass(frozen=True)
class RiskScenario:
    pga_level: str = "10%"
    fault_lambda_km: float = 200.0
    amplification_method: str = "boore"
    analysis_type: str = "base_case"

    @property
    def pga_label(self) -> str:
        return self.pga_level if self.pga_level.endswith("%") else f"{self.pga_level}%"


def _normalize_age(age: object) -> str:
    if pd.isna(age):
        return ""
    return str(age).strip().lower()


def age_weight(age: object) -> float:
    a = _normalize_age(age)
    if not a:
        return 1.0
    if a.startswith("qh"):
        return 4.0
    if a in {"q", "qp3"} or "qp3" in a:
        return 3.0
    if "qp" in a:
        return 2.0
    return 1.0


def fea_weight(fea_en: object, fea_ch: object = None) -> float:
    text = f"{fea_en or ''} {fea_ch or ''}".lower()
    if any(k in text for k in ("buried", "inferred", "隐", "推断")):
        return 4.0
    if any(k in text for k in ("oblique", "complex", "斜", "走滑正断", "走滑逆断")):
        return 3.0
    if any(
        k in text
        for k in (
            "normal", "reverse", "strike", "lateral", "left", "right",
            "正断", "逆断", "走滑", "左旋", "右旋",
        )
    ):
        return 2.0
    return 1.0


def fais(age: object, fea_en: object, fea_ch: object = None) -> float:
    return age_weight(age) * 1.5 + fea_weight(fea_en, fea_ch)


def _acute_alpha_deg(site_azimuth_rad: float, fault_azimuth_rad: float) -> float:
    diff = abs(site_azimuth_rad - fault_azimuth_rad)
    alpha = math.degrees(diff) % 180.0
    if alpha > 90.0:
        alpha = 180.0 - alpha
    return alpha


def _segment_strike_azimuth(line: LineString, near_pt: Point) -> float:
    coords = list(line.coords)
    if len(coords) < 2:
        return 0.0
    distances = [Point(c).distance(near_pt) for c in coords]
    idx = int(np.argmin(distances))
    if idx == len(coords) - 1:
        x0, y0 = coords[idx - 1]
        x1, y1 = coords[idx]
    else:
        x0, y0 = coords[idx]
        x1, y1 = coords[idx + 1]
    return math.atan2(y1 - y0, x1 - x0)


def compute_alpha_deg(site_proj: Point, fault_geom_proj: LineString) -> float:
    site_pt, fault_pt = nearest_points(site_proj, fault_geom_proj)
    site_az = math.atan2(fault_pt.y - site_pt.y, fault_pt.x - site_pt.x)
    fault_az = _segment_strike_azimuth(fault_geom_proj, fault_pt)
    return _acute_alpha_deg(site_az, fault_az)


def sfp_score_array(
    fais_vals: np.ndarray,
    distance_km: np.ndarray,
    alpha_deg: np.ndarray,
    lambda_km: float,
) -> np.ndarray:
    directional = 0.5 * (1.0 + np.cos(np.radians(alpha_deg)))
    raw = fais_vals * np.exp(-distance_km / lambda_km) * directional
    return raw / FAIS_MAX * 10.0


def spga_score(pga_value: np.ndarray | float, pga_max: float) -> np.ndarray:
    if pga_max <= 0:
        return np.zeros_like(np.asarray(pga_value, dtype=np.float64))
    return np.asarray(pga_value, dtype=np.float64) / pga_max * 10.0


def _f_from_vs30(vs30: np.ndarray) -> np.ndarray:
    vs30 = np.maximum(vs30.astype(np.float64), 1.0)
    return np.clip((vs30 / BOORE_VREF) ** BOORE_B, 0.5, 3.0)


def ssa_from_vs30(vs30: np.ndarray, method: str) -> np.ndarray:
    method = method.lower()
    if method == "none":
        f = np.full_like(vs30, 1.0, dtype=np.float64)
    elif method in {"gb50011", "gb_50011", "code"}:
        f = np.ones_like(vs30, dtype=np.float64)
        f = np.where(vs30 >= 500.0, 0.8, f)
        f = np.where((vs30 >= 250.0) & (vs30 < 500.0), 0.9, f)
        f = np.where((vs30 >= 150.0) & (vs30 < 250.0), 1.0, f)
        f = np.where(vs30 < 150.0, 1.2, f)
    else:
        f = _f_from_vs30(vs30)
    return np.clip((f - 0.5) / 2.5 * 10.0, 0.0, 10.0)


def final_risk_score(
    spga: np.ndarray,
    sfp: np.ndarray,
    ssa: np.ndarray,
    w1: float = W1,
    w2: float = W2,
    w3: float = W3,
) -> np.ndarray:
    return w1 * spga + w2 * sfp + w3 * ssa


# Weight sensitivity (Risk Assessment Redo.docx, last table)
# Columns: w1 (SPGA), w2 (SFP), w3 (SSA)
WEIGHT_SCENARIOS: list[tuple[str, str, tuple[float, float, float]]] = [
    ("sensitivity_weight_tectonic_risk.csv", "Tectonic Dominance", (0.30, 0.50, 0.20)),
    ("sensitivity_weight_geotechnical_risk.csv", "Geotechnical Dominance", (0.40, 0.20, 0.40)),
    ("sensitivity_weight_shaking_risk.csv", "Absolute Shaking", (0.60, 0.30, 0.10)),
]


def reweight_risk_table(
    base: pd.DataFrame,
    w1: float,
    w2: float,
    w3: float,
    analysis_type: str = "weight_sensitivity",
) -> pd.DataFrame:
    """Recompute final_risk_score from stored SPGA/SFP/SSA component columns."""
    out = base.copy()
    spga = out["pga_score"].to_numpy(dtype=np.float64)
    sfp = out["fault_score"].to_numpy(dtype=np.float64)
    ssa = out["site_vulnerability"].to_numpy(dtype=np.float64)
    final_r = final_risk_score(spga, sfp, ssa, w1=w1, w2=w2, w3=w3)
    out["final_risk_score"] = final_r
    out["risk_category"] = risk_category(final_r)
    out["analysis_type"] = analysis_type
    out["w1"] = w1
    out["w2"] = w2
    out["w3"] = w3
    return out


def srf_listed_colormap(vmin: float = 0, vmax: float = 10):
    """
    Discrete SRF colormap aligned with risk_category bands:
      Low (0–3): greens → Moderate (3–6): blues → High / Very High (>6): reds.
    Returns (cmap, norm, boundaries).
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    vmax = max(float(vmax), 6.0)
    boundaries = np.arange(vmin, vmax + 2)
    colors = []
    for val in boundaries[:-1]:
        if val < 3:
            frac = val / 3 if 3 > 0 else 0
            colors.append(plt.cm.Greens(0.35 + frac * 0.55))
        elif val < 6:
            frac = (val - 3) / 3
            colors.append(plt.cm.Blues(0.40 + frac * 0.50))
        else:
            frac = (val - 6) / max(1.0, vmax - 6)
            colors.append(plt.cm.Reds(0.35 + frac * 0.60))
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries, cmap.N)
    return cmap, norm, boundaries


# Display name after manuscript revision (was "Seismic risk factor (SRF)")
SSI_COLORBAR_LABEL = "Seismic Screening Index (SSI)"
SSI_RISK_LEGEND = (
    "0≤SSI≤3 (Low Risk) | "
    "3<SSI≤6 (Moderate Risk) | "
    "6<SSI≤9 (High Risk) | "
    "SSI>9 (Very High Risk)"
)


def style_srf_colorbar(
    fig,
    cbar,
    label: str | None = None,
    fontsize_label: float | None = None,
    fontsize_tick: float = 7.0,
    label_dx: float = 0.04,
    *,
    add_label: bool = True,
):
    """Style a vertical SSI colorbar; optionally place the axis label past the ticks."""
    cbar.ax.tick_params(labelsize=fontsize_tick, pad=4)
    cbar.outline.set_linewidth(0.6)
    try:
        cbar.dividers.set_linewidth(0.4)
    except Exception:
        pass
    if add_label:
        if label is None:
            label = SSI_COLORBAR_LABEL
        if fontsize_label is None:
            fontsize_label = fontsize_tick
        pos = cbar.ax.get_position()
        fig.text(
            pos.x1 + label_dx,
            0.5 * (pos.y0 + pos.y1),
            label,
            rotation=90,
            va="center",
            ha="left",
            fontsize=fontsize_label,
        )
    return cbar


def add_shared_ssi_colorbar_label(
    fig,
    cbars: list,
    fontsize_label: float,
    label_dx: float = 0.04,
    label: str = SSI_COLORBAR_LABEL,
):
    """One vertical SSI label centered across multiple row colorbars."""
    if not cbars:
        return
    positions = [cbar.ax.get_position() for cbar in cbars]
    x = max(p.x1 for p in positions) + label_dx
    y0 = min(p.y0 for p in positions)
    y1 = max(p.y1 for p in positions)
    fig.text(
        x,
        0.5 * (y0 + y1),
        label,
        rotation=90,
        va="center",
        ha="left",
        fontsize=fontsize_label,
    )


def risk_category(values: np.ndarray) -> np.ndarray:
    out = np.full(len(values), "Low", dtype=object)
    out[(values > 3) & (values <= 6)] = "Moderate"
    out[(values > 6) & (values <= 9)] = "High"
    out[values > 9] = "Very High"
    return out


def load_storage_sites() -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for path, stype in ((DSA_PIPELINE, "DSA"), (EOR_PIPELINE, "EOR")):
        df = pd.read_excel(path, usecols=["X", "Y", "Storage Potential"])
        df = df.rename(columns={"X": "lon", "Y": "lat", "Storage Potential": "storage_potential"})
        df["storage_type"] = stype
        parts.append(df)
    out = pd.concat(parts, ignore_index=True)
    return out.dropna(subset=["lon", "lat", "storage_potential"])


@lru_cache(maxsize=1)
def _load_fault_segments(cafd_path: str) -> gpd.GeoDataFrame:
    usecols = ["id", "X (lon)", "Y (lat)", "vertex_index", "AGE", "Fea_En", "Fea_Ch"]
    raw = pd.read_csv(cafd_path, usecols=usecols)
    raw = raw.dropna(subset=["id", "X (lon)", "Y (lat)"])
    segments: list[dict] = []
    for seg_id, grp in raw.groupby("id", sort=False):
        grp = grp.sort_values("vertex_index")
        if len(grp) < 2:
            continue
        line = LineString(list(zip(grp["X (lon)"], grp["Y (lat)"])))
        if line.length == 0:
            continue
        row = grp.iloc[0]
        segments.append(
            {
                "seg_id": seg_id,
                "geometry": line,
                "AGE": row["AGE"],
                "Fea_En": row["Fea_En"],
                "Fea_Ch": row["Fea_Ch"],
                "fais": fais(row["AGE"], row["Fea_En"], row["Fea_Ch"]),
            }
        )
    gdf = gpd.GeoDataFrame(segments, geometry="geometry", crs="EPSG:4326")
    return gdf.to_crs(CHINA_ALBERS)


def _nearest_fault_metrics(
    lons: np.ndarray,
    lats: np.ndarray,
    faults: gpd.GeoDataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sites = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(lons, lats),
        crs="EPSG:4326",
    ).to_crs(CHINA_ALBERS)

    tree = STRtree(faults.geometry.values)
    distances_km = np.zeros(len(sites), dtype=np.float64)
    alphas = np.zeros(len(sites), dtype=np.float64)
    fais_vals = np.zeros(len(sites), dtype=np.float64)

    geoms = faults.geometry.values
    fais_arr = faults["fais"].to_numpy(dtype=np.float64)

    for i, pt in enumerate(sites.geometry):
        idx = tree.nearest(pt)
        fault_geom = geoms[idx]
        distances_km[i] = pt.distance(fault_geom) / 1000.0
        alphas[i] = compute_alpha_deg(pt, fault_geom)
        fais_vals[i] = fais_arr[idx]

    return distances_km, alphas, fais_vals


@lru_cache(maxsize=8)
def _load_grid_values(csv_path: str, value_col: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    lon_col = "Long" if "Long" in df.columns else "X"
    lat_col = "Lat" if "Lat" in df.columns else "Y"
    if value_col not in df.columns:
        matches = [c for c in df.columns if value_col.replace("%", "") in c.replace("%", "")]
        if not matches:
            raise ValueError(f"Column {value_col!r} not found in {csv_path}")
        value_col = matches[0]
    lons = df[lon_col].to_numpy(dtype=np.float64)
    lats = df[lat_col].to_numpy(dtype=np.float64)
    vals = df[value_col].to_numpy(dtype=np.float64)
    return lons, lats, vals


def _sample_grid(lons: np.ndarray, lats: np.ndarray, grid_lon, grid_lat, grid_val) -> np.ndarray:
    tree = BallTree(np.radians(np.column_stack([grid_lat, grid_lon])), metric="haversine")
    q = np.radians(np.column_stack([lats, lons]))
    _, idx = tree.query(q, k=1)
    return grid_val[idx[:, 0]]


def compute_risk_table(scenario: RiskScenario) -> pd.DataFrame:
    sites = load_storage_sites()
    faults = _load_fault_segments(str(CAFD_FAULTS_CSV))

    dist_km, alpha, fais_vals = _nearest_fault_metrics(
        sites["lon"].to_numpy(),
        sites["lat"].to_numpy(),
        faults,
    )
    sfp = sfp_score_array(fais_vals, dist_km, alpha, scenario.fault_lambda_km)

    pga_path = PGA_CSV_BY_LEVEL[scenario.pga_label]
    pga_col = f"PE50a{scenario.pga_label.replace('%', '')}%"
    g_lon, g_lat, g_pga = _load_grid_values(str(pga_path), pga_col)
    pga_site = _sample_grid(
        sites["lon"].to_numpy(), sites["lat"].to_numpy(), g_lon, g_lat, g_pga,
    )
    spga = spga_score(pga_site, float(np.max(g_pga)))

    v_lon, v_lat, v_vs30 = _load_grid_values(str(VS30_CSV), "Vs30")
    vs30_site = _sample_grid(
        sites["lon"].to_numpy(), sites["lat"].to_numpy(), v_lon, v_lat, v_vs30,
    )
    ssa = ssa_from_vs30(vs30_site, scenario.amplification_method)
    final_r = final_risk_score(spga, sfp, ssa)

    out = sites.copy()
    out["pga_score"] = spga
    out["fault_score"] = sfp
    out["site_vulnerability"] = ssa
    out["final_risk_score"] = final_r
    out["pga_level"] = scenario.pga_label
    out["fault_lambda"] = scenario.fault_lambda_km
    out["amplification_method"] = scenario.amplification_method
    out["analysis_type"] = scenario.analysis_type
    out["risk_category"] = risk_category(final_r)
    out["nearest_fault_km"] = dist_km
    out["fault_alpha_deg"] = alpha
    out["fais"] = fais_vals
    return out


DEFAULT_SCENARIOS: list[tuple[str, RiskScenario]] = [
    ("final_seismic_risk_factor_base_case.csv", RiskScenario()),
    ("sensitivity_pga_0.5_risk.csv", RiskScenario(pga_level="0.5%", analysis_type="pga_sensitivity")),
    ("sensitivity_pga_1_risk.csv", RiskScenario(pga_level="1%", analysis_type="pga_sensitivity")),
    ("sensitivity_pga_2_risk.csv", RiskScenario(pga_level="2%", analysis_type="pga_sensitivity")),
    ("sensitivity_pga_5_risk.csv", RiskScenario(pga_level="5%", analysis_type="pga_sensitivity")),
    ("sensitivity_pga_63_risk.csv", RiskScenario(pga_level="63%", analysis_type="pga_sensitivity")),
    ("sensitivity_fault_lambda_50_risk.csv", RiskScenario(fault_lambda_km=50, analysis_type="lambda_sensitivity")),
    ("sensitivity_fault_lambda_100_risk.csv", RiskScenario(fault_lambda_km=100, analysis_type="lambda_sensitivity")),
    ("sensitivity_fault_lambda_300_risk.csv", RiskScenario(fault_lambda_km=300, analysis_type="lambda_sensitivity")),
    ("sensitivity_fault_lambda_500_risk.csv", RiskScenario(fault_lambda_km=500, analysis_type="lambda_sensitivity")),
    ("sensitivity_amplification_gb50011_risk.csv", RiskScenario(amplification_method="gb50011", analysis_type="amplification_sensitivity")),
    ("sensitivity_amplification_none_risk.csv", RiskScenario(amplification_method="none", analysis_type="amplification_sensitivity")),
]
