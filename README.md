# CCS Figures

Python tooling for Carbon Capture and Storage (CCS) figures and supporting tables for China-focused analysis (cost maps, seismic risk, source–sink matching, Sankey flows, blue hydrogen, and efficiency radials).

## Setup

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate

# Linux / macOS
source .venv/bin/activate

pip install -r requirements.txt
```

Python 3.11+ recommended.

## Repository layout

```
CCS Figures/
├── charts/              # Figure outputs only (PNG, HTML)
├── data/                # All inputs and tabular outputs (CSV, XLSX, shapefiles)
├── sankey/              # Optional standalone Sankey web assets
├── config.py            # Shared paths, fonts, Plotly PNG DPI
├── run_all_figures.py   # Batch runner
├── requirements.txt
└── figure*.py           # Figure scripts (+ helper modules)
```

### `charts/` — figures only

PNG maps/plots and Plotly HTML Sankeys (`DSA_Sankey_Plot.html`, `EOR_Sankey_Plot.html`).  
Use `chart_path("name.png")` from `config.py`.

### `data/` — inputs and tables

| Path | Contents |
|------|----------|
| `DSA-Pipeline.xlsx`, `DSA-Tank.xlsx`, `EOR_Pipeline.xlsx`, `EOR_Tank.xlsx` | Storage / cost grids |
| `industries.csv`, `china_gas_fields.csv` | Facilities and gas fields (`industries.csv` is on Zenodo; see below) |
| `Carbon_DSA.xlsx`, `Carbon_EOR.xlsx`, `Emission_Storage.xlsx` | Carbon / emission series |
| `natural_earth/` | Admin0 / Admin1 boundaries |
| `Risk_Assessment/` | Seismic risk CSVs and weight-sensitivity tables (the large CAFD400 geom CSV is on Zenodo) |
| `ss_scores/` | Sink suitability score CSVs (on Zenodo) |
| `Ex-En-Subsystem.xlsx`, `Generation streams-System *.xlsx` | Efficiency radial inputs |
| `*.csv` at data root | Derived tables (Sankey flows, supplementary cost stats, blue-H₂ assignments, industry stats, …) |

Use `data_path("name.csv")` for tabular outputs. Do not park CSVs under `charts/`.

Large inputs are omitted from git. Download [ccs_large_files.zip](https://doi.org/10.5281/zenodo.22110783) from Zenodo and extract so these paths exist:

- `data/industries.csv`
- `data/ss_scores/`
- `data/Risk_Assessment/CAFD400_V2023_1-Reprojected-5km-Intersection-Aggregate-GPSG-4326-Geom.csv`

**DOI:** [10.5281/zenodo.22110783](https://doi.org/10.5281/zenodo.22110783)

## Running figures

```bash
# Single figure
python figure3a.py

# All scripts listed in run_all_figures.py
python run_all_figures.py
```

Outputs land in `charts/`; derived CSVs land in `data/` (or `data/Risk_Assessment/` for risk tables).

## Main modules

| Module | Role |
|--------|------|
| `config.py` | Paths (`DATA_DIR`, `CHARTS_DIR`), manuscript fonts, Plotly DPI |
| `cost_maps.py` / `cost_diff_maps.py` | Cost and cost-difference maps |
| `carbon_contributions.py` | Carbon contribution bars |
| `seismic_risk.py` / `compute_seismic_risk.py` | Seismic risk scoring and recompute |
| `blue_hydrogen_smr.py` / `figure6_blue_hydrogen.py` | Dual-gate SMR map and cost-supply curve |
| `didpa.py` | Source–sink matching / SS scores |
| `plotly_export.py` | Plotly HTML export with PNG DPI metadata |
| `province_style.py` | Shared province styling |
| `export_sankey_data.py` | Province flow tables for Sankey |

## Notes

- Grid costs are in **USD/tCO₂**.
- Manuscript text sizing uses `manuscript_font_bundle(fig.get_figwidth())` so labels read ~7 pt at a 6.5″ figure width.
- Spec Word docs at the repo root (`Blue hydrogen redo.docx`, `Risk Assessment Redo.docx`, …) document methods; they are not required to run plots if the Excel/CSV inputs are present.
- Large data files listed in `.gitignore` are on Zenodo: [10.5281/zenodo.22110783](https://doi.org/10.5281/zenodo.22110783).
