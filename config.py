"""Shared paths for CCS Figures scripts."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"

# Grid-level DSA / EOR cost and storage files (USD/tCO₂)
DSA_PIPELINE = DATA_DIR / "DSA-Pipeline.xlsx"
DSA_TANK = DATA_DIR / "DSA-Tank.xlsx"
EOR_PIPELINE = DATA_DIR / "EOR_Pipeline.xlsx"
EOR_TANK = DATA_DIR / "EOR_Tank.xlsx"

# Grid storage potential + injectivity (X, Y; Storage Potential / Injection Rate; R_yyyy, S_yyyy)
STORAGE_POTENTIAL_DSA_PIPELINE = DSA_PIPELINE
STORAGE_POTENTIAL_DSA_TANK = DSA_TANK
STORAGE_POTENTIAL_EOR_PIPELINE = EOR_PIPELINE
STORAGE_POTENTIAL_EOR_TANK = EOR_TANK
# Backward-compatible aliases (pipeline-only)
STORAGE_POTENTIAL_DSA = DSA_PIPELINE
STORAGE_POTENTIAL_EOR = EOR_PIPELINE

# All grid sink workbooks for SS / matching (pipeline + tank, DSA + EOR)
STORAGE_POTENTIAL_GRIDS = (
    STORAGE_POTENTIAL_DSA_PIPELINE,
    STORAGE_POTENTIAL_DSA_TANK,
    STORAGE_POTENTIAL_EOR_PIPELINE,
    STORAGE_POTENTIAL_EOR_TANK,
)

# Supporting inputs that remain in data/
INDUSTRIES_CSV = DATA_DIR / "industries.csv"
CARBON_DSA = DATA_DIR / "Carbon_DSA.xlsx"
CARBON_EOR = DATA_DIR / "Carbon_EOR.xlsx"
EMISSION_STORAGE = DATA_DIR / "Emission_Storage.xlsx"
NE_ADMIN0 = DATA_DIR / "natural_earth" / "admin0"
NE_ADMIN1 = DATA_DIR / "natural_earth" / "admin1"
CHINA_GAS_FIELDS = DATA_DIR / "china_gas_fields.csv"

# Blue hydrogen Dual-Gate SMR (Blue hydrogen redo.docx; no 1.3 winding factor on distances)
BLUE_H2_INDUSTRIES = ("MC19", "PM19")
BLUE_H2_RADIUS_KM = 250
BLUE_H2_RADIUS_M = BLUE_H2_RADIUS_KM * 1000
# Pipeline distance cost weights (Eq. D terms); normalize CO2 leg to 1.0
W_CO2 = 1.0
W_H2 = 1.0
W_GAS = 1.0
# H_avg = S_H2 / (L * r); cost-supply x-axis capped at max H_2050 (Table 4)
BLUE_H2_LIFESPAN_YEARS = 30
BLUE_H2_CO2_TO_H2_RATIO = 10.0
BLUE_H2_H2050_XMAX = 7.6

COST_YEARS = [2025, 2030, 2035, 2040, 2045, 2050]
COST_SCENARIOS = ("min", "avg", "max")
CNY_TO_USD = 1 / 7.183

# Grid costs are in USD/tCO2.
DATA_COST_UNIT = "USD/tCO₂"
# Backward-compatible aliases
DATA2_COST_UNIT = DATA_COST_UNIT
LEGACY_COST_UNIT = "USD/km"

CHARTS_DIR = ROOT / "charts"
SS_SCORES_DIR = DATA_DIR / "ss_scores"
RISK_DIR = DATA_DIR / "Risk_Assessment"

# Manuscript text rule: when a figure is placed at MANUSCRIPT_FIG_WIDTH_IN,
# on-figure text should appear as ~MANUSCRIPT_TEXT_PT (journal-style 7 pt).
MANUSCRIPT_TEXT_PT = 7.0
MANUSCRIPT_FIG_WIDTH_IN = 6.5  # typical double-column / full-text width in Word


def manuscript_fontsize(design_width_in: float, target_pt: float = MANUSCRIPT_TEXT_PT) -> float:
    """Design fontsize so text reads as target_pt at manuscript width.

    print_pt ≈ design_pt * (manuscript_width / design_width)
    ⇒ design_pt = target_pt * (design_width / manuscript_width)
    """
    if design_width_in <= 0:
        return float(target_pt)
    return float(target_pt) * (float(design_width_in) / MANUSCRIPT_FIG_WIDTH_IN)


def manuscript_font_bundle(design_width_in: float, target_pt: float = MANUSCRIPT_TEXT_PT) -> dict:
    """Font sizes for a figure of the given design width (inches).

    Prefer this over module-level FONT_* when figsize is known — especially for
    wide multi-panel canvases (e.g. 20 in), where a fixed 10 in calibration
    would print too small at manuscript width.
    """
    body = round(manuscript_fontsize(design_width_in, target_pt), 1)
    return {
        "body": body,
        "label": body,
        "tick": body,
        "title": body,
        "legend": body,
        "legend_marker": max(6.0, round(body - 2.0, 1)),
        "colorbar_label": body,
        "colorbar_tick": body,
    }


# Matplotlib defaults assume ~10 in wide China-map canvases (figsize ≈ 10).
# Prefer manuscript_font_bundle(fig.get_figwidth()) inside plot functions.
_REF_MAP_WIDTH_IN = 10.0
_MF = manuscript_font_bundle(_REF_MAP_WIDTH_IN)
FONT_AXIS_LABEL = _MF["label"]
FONT_AXIS_TICK = _MF["tick"]
FONT_COLORBAR_LABEL = _MF["colorbar_label"]
FONT_COLORBAR_TICK = _MF["colorbar_tick"]
FONT_PANEL_TITLE = _MF["title"]
FONT_PANEL_LABEL = FONT_PANEL_TITLE  # alias used by cost map scripts
FONT_LEGEND = _MF["legend"]
FONT_LEGEND_MARKER = _MF["legend_marker"]

# Plotly HTML toolbar PNG export
PLOTLY_PNG_DPI = 300
PLOTLY_PNG_SCALE = PLOTLY_PNG_DPI / 96

# Seismic risk inputs (Risk Assessment Redo / EE-ART manuscript Eq. 5–11)
CAFD_FAULTS_CSV = RISK_DIR / "CAFD400_V2023_1-Reprojected-5km-Intersection-Aggregate-GPSG-4326-Geom.csv"
VS30_CSV = RISK_DIR / "VS30_China_5km_Interpolated.csv"
PGA_CSV_BY_LEVEL = {
    "0.5%": RISK_DIR / "PGA-PE50a0.5-IDW-ESRI-102012-5km-2.csv",
    "1%": RISK_DIR / "PGA-PE50a1-IDW-ESRI-102012-5km.csv",
    "2%": RISK_DIR / "PGA-PE50a2-IDW-ESRI-102012-5km.csv",
    "5%": RISK_DIR / "PGA-PE50a5-IDW-ESRI-102012-5km.csv",
    "10%": RISK_DIR / "PGA-PE50a10-IDW-ESRI-102012-5km.csv",
    "63%": RISK_DIR / "PGA-PE50a63-IDW-ESRI-102012-5km.csv",
}


def _disable_matplotlib_grids() -> None:
    """Publication default: no axis gridlines on any matplotlib figure."""
    try:
        import matplotlib as mpl

        mpl.rcParams["axes.grid"] = False
    except ImportError:
        pass


_disable_matplotlib_grids()


def ensure_charts_dir() -> Path:
    CHARTS_DIR.mkdir(parents=True, exist_ok=True)
    return CHARTS_DIR


def ensure_data_dir() -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    return DATA_DIR


def chart_path(filename: str) -> str:
    """Return absolute path for a figure output inside charts/."""
    return str(ensure_charts_dir() / filename)


def data_path(filename: str) -> str:
    """Return absolute path for a tabular/data output inside data/."""
    return str(ensure_data_dir() / filename)
