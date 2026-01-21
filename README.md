# China-CCS-Project
# CCS Figures

Python scripts for generating figures related to Carbon Capture and Storage (CCS) analysis in China.

## Setup

### Prerequisites
- Python 3.11 or higher
- Virtual environment recommended

### Installation

```bash
# Create virtual environment (from project root)
python -m venv .venv

# Activate virtual environment
# Windows:
.venv\Scripts\activate
# Linux/Mac:
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## Data Directory Structure

All scripts expect data files in the `./data/` directory:

```
CCS Figures/
├── data/
│   ├── Carbon_DSA.xlsx          # DSA carbon contribution data
│   ├── Carbon_EOR.xlsx          # EOR carbon contribution data
│   ├── DSA-Pipeline.xlsx        # DSA pipeline storage data
│   ├── DSA-Tank.xlsx            # DSA tank storage data
│   ├── EOR_Pipeline.xlsx        # EOR pipeline storage data
│   ├── EOR_Tank.xlsx            # EOR tank storage data
│   ├── Emission_Storage.xlsx    # Emission and storage projections
│   ├── china_gas_fields.csv     # Gas field locations
│   ├── industries.csv           # Industrial facility locations
│   ├── natural_earth/           # Shapefile directories
│   │   ├── admin0/              # Country boundaries
│   │   └── admin1/              # Province boundaries
│   └── Risk_Assessment/         # Seismic risk data files
│       ├── CAFD400-SRF-5km-2.csv
│       ├── final_seismic_risk_factor_base_case.csv
│       ├── PGA-*.csv            # Peak Ground Acceleration files
│       ├── sensitivity_*.csv    # Sensitivity analysis files
│       └── VS30_China_5km_Interpolated.csv
```

## Scripts Overview

### Core Modules

| Module | Description |
|--------|-------------|
| `cost_maps.py` | Cost map plotting functions |
| `carbon_contributions.py` | Carbon contribution bar chart functions |

### Main Figures

| Script | Description | Output |
|--------|-------------|--------|
| `figure1.py` | Injection maps (Storage Potential & Injection Rate) | `injection_maps.png` |
| `figure2.py` | Cost maps (Tank, avg) | `cost_maps_tank_avg.png` |
| `figure3a.py` | Seismic risk factor map | `seismic_risk_factor_map.png` |
| `figure3b.py` | Risk vs storage bubble chart | `risk_storage_bubble.png` |
| `figure4a.py` | Storage potential and industries map | `storage_industry_map.png` |
| `figure4b.py` | Storage vs injection scatter plot | `storage_vs_injection.png` |
| `figure4c.py` | Storage comparison (DSA vs EOR) | `storage_comparison.png` |
| `figure4d.py` | Injection vs unabated emissions | `injection_unabated.png` |
| `figure5a-b.py` | CO2 flow Sankey diagrams | `DSA_Sankey_Plot.html`, `EOR_Sankey_Plot.html` |
| `figure6.py` | Storage potential and gas fields map | `storage_gas_map.png` |
| `figure7.py` | Carbon bubble charts | `Carbon_DSA_bubble_chart.png`, `Carbon_EOR_bubble_chart.png` |

### Supplementary Figures

| Script | Description | Output |
|--------|-------------|--------|
| `figure1_sup.py` | Storage and industries (supplementary) | `storage_industry_map_sup.png` |
| `figure2-3_sup.py` | Cost maps (Pipeline, min) | `cost_maps_pipeline_min.png` |
| `figure4-5_sup.py` | Cost maps (Tank, min) | `cost_maps_tank_min.png` |
| `figure6-7_sup.py` | Cost maps (Pipeline, avg) | `cost_maps_pipeline_avg.png` |
| `figure8-9_sup.py` | Cost maps (Tank, avg) | `cost_maps_tank_avg_sup.png` |
| `figure10-11_sup.py` | Cost maps (Pipeline, max) | `cost_maps_pipeline_max.png` |
| `figure12-13_sup.py` | Cost maps (Tank, max) | `cost_maps_tank_max.png` |
| `figure14_sup.py` | Seismic risk score (SRS) map | `seismic_risk_srs_map.png` |
| `figure15_sup.py` | Historical seismicity map | `historical_seismicity_map.png` |
| `figure16_sup.py` | Multi-panel PGA maps | `pga_multi_panel.png` |
| `figure17_sup.py` | Sensitivity analysis - PGA thresholds | `sensitivity_pga.png` |
| `figure18_sup.py` | Sensitivity analysis - Site amplification | `sensitivity_amplification.png` |
| `figure19_sup.py` | Sensitivity analysis - Fault lambda | `sensitivity_fault_lambda.png` |
| `figure20_sup.py` | Storage and industries (multi-plot) | `storage_industry_*.png` |
| `figure21a_sup.py` | Carbon contributions (DSA) | `carbon_contributions_dsa.png` |
| `figure21b_sup.py` | Carbon contributions (EOR) | `carbon_contributions_eor.png` |
| `figure22_sup.py` | VS30 map | `vs30_map.png` |

## Running Scripts

Each script can be run independently from the `CCS Figures` directory:

```bash
cd "CCS Figures"
python figure1.py
python figure2.py
# etc.
```

Or from the project root:

```bash
python "CCS Figures/figure1.py"
```

## Map Projections

Most geospatial figures use:
- **Input CRS**: EPSG:4326 (WGS84)
- **Target CRS**: ESRI:102012 (Asia Lambert Conformal Conic)

## Color Conventions

- **Green**: DSA (Deep Saline Aquifer) storage
- **Orange/Gray**: EOR (Enhanced Oil Recovery) storage
- **Red**: Negative or clipped values
- **Grey**: Special zero conditions
- **Blue hatched**: Overlap/buffer zones
- **Viridis/Turbo colormaps**: Continuous value ranges
