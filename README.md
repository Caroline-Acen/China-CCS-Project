# China CCS figures

Python scripts that reproduce the manuscript and supplementary figures for the China CCS analysis (cost maps, seismic screening, source–sink matching, Sankey flows, blue hydrogen, and cost comparison).

## Setup

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate

# Linux / macOS
source .venv/bin/activate

pip install -r requirements.txt
```

Python 3.11+ is recommended. Several maps need a working GeoPandas stack (GEOS/GDAL).

## Large data (Zenodo)

These paths are gitignored because of size. Download [ccs_large_files.zip](https://doi.org/10.5281/zenodo.22110783) and extract so they exist:

- `data/industries.csv`
- `data/ss_scores/`
- `data/Risk_Assessment/CAFD400_V2023_1-Reprojected-5km-Intersection-Aggregate-GPSG-4326-Geom.csv`

**DOI:** [10.5281/zenodo.22110783](https://doi.org/10.5281/zenodo.22110783)

Without them, cost maps, sink-suitability Figure 5b, industry maps, and the fault map (Supplementary Figure 2) will fail. Dual-Gate assignment CSVs are already in `data/`; Figure 9 reuses them unless you pass `--recompute-assignments`.

### Source–sink distances

Industry grids are matched to storage cells within **250 km** (great-circle). For each grid the files store the distance to the sink that maximises the sink-suitability score.

| Path | Contents |
|---|---|
| `data/sink_suitability_summary.csv` | Per year (2025–2050): count of grids and `best_distance_km_mean`, `_median`, `_p85`, `_max` |
| `data/source_sink_model_haul_km.csv` | The same yearly stats, plus pooled all-years mean and median |
| `data/ss_scores/` (Zenodo) | Per-category and combined CSVs with column `Best_Distance_km` |
| `data/gasfield_gate_smr_assignments.csv` | Figure 9 Dual-Gate assignments: `d_co2_km`, `d_h2_km` |
| `data/sink_gate_smr_assignments.csv` | Figure 9 Dual-Gate assignments: `d_gas_km`, `d_h2_km` |

Supplementary Figure 19 (historical seismicity) also needs `data/historical_events.xlsx` (latitude, longitude, magnitude). That file is not in the repository or the Zenodo zip; the script skips if it is absent.

## Reproduce all figures

```bash
python run_all_figures.py
```

That:

1. Writes native outputs into `charts/` (descriptive filenames).
2. Copies/upscales manuscript names into `charts/high_resolution/` (`figure1.png`, `sup_figure1.png`, …).
3. Overwrites `high_resolution/figure4b.png` (3150×3000) and `figure9b.png` (3000×3000) at the journal pixel sizes.

Single figure:

```bash
python figure4a.py
python export_high_resolution_charts.py
```

## Output layout

```
charts/                 # native PNG/HTML from each script
charts/high_resolution/ # paper IDs, ≥3000 px on each edge (except SF24, copied as-is)
data/                   # inputs and derived tables
```

`charts/` keeps descriptive PNG/HTML names. `charts/high_resolution/` uses manuscript IDs. Python scripts use the same IDs (`python figure4b.py`, `python sup_figure1.py`). Where two panels share a script, a short wrapper exists for the other ID (e.g. `figure8.py` runs `figure7.py`).

| Paper file | Source in `charts/` |
|---|---|
| `figure1.png` | `EOR+DSA pipeline avg.png` |
| `figure2.png` | `EOR+DSA tank avg.png` |
| `figure3.png` | `Cost difference (EOR+DSA tank avg - EOR+DSA pipeline avg).png` |
| `figure4a.png` | `seismic_risk_factor_map.png` |
| `figure4b.png` | `risk_storage_bubble.png` (journal pixels in high-res) |
| `figure5a.png` | `storage_industry_map.png` |
| `figure5b.png` | `sink_suitability_map_panels.png` |
| `figure6a.png` / `figure6b.png` / `figure6c.png` | `storage_comparison_dsa.png`, `storage_comparison_eor.png`, `injection_unabated.png` |
| `figure7.html` / `figure8.html` | `DSA_Sankey_Plot.html`, `EOR_Sankey_Plot.html` |
| `figure9a.png` / `figure9b.png` | `gasfield_sink_combined_gate_map.png`, `hydrogen_hub_cost_supply_curve.png` |
| `figure10a.png` / `figure10b.png` | `Carbon_DSA_bubble_chart.png`, `Carbon_EOR_bubble_chart.png` |
| `figure11.png` | `figure11.png` (250 km radius). Mean/median hauls: `figure11_mean.png`, `figure11_median.png` |
| `sup_figure1.png` | `pga_multi_panel.png` |
| `sup_figure2.png` | `seismic_faults_by_name_map_v2.png` |
| `sup_figure3.png` | `vs30_map.png` |
| `sup_figure4.png` | copy of `figure5a.png` |
| `sup_figure5.png` | `storage_gas_map.png` |
| `sup_figure6.png`–`sup_figure17.png` | DSA/EOR pipeline and tank min/avg/max cost maps |
| `sup_figure18.png` | `injection_maps.png` |
| `sup_figure19.png` | `historical_seismicity_map.png` (optional data) |
| `sup_figure20.png`–`sup_figure23.png` | PGA / amplification / fault-λ / weight sensitivity |
| `sup_figure24a.png` / `sup_figure24b.png` | `sup_figure24a.png`, `sup_figure24b.png` (copied, not upscaled) |

## Figure 11 units

Supplementary Table 15: EU values are converted with **1 USD = 0.87 EUR**. USA and EU transport are already per-ton totals. China transport is 0.11–0.19 USD t⁻¹ km⁻¹; `figure11.py` writes `figure11_mean.png` and `figure11_median.png` using the pooled source–sink distances in `data/source_sink_model_haul_km.csv`. `figure11.png` still uses the 250 km matching radius.

## Notes

- Grid costs are **USD/tCO₂**.
- On-figure type uses `manuscript_font_bundle(...)` so labels read ~7 pt at a 6.5″ manuscript width.
- High-resolution carbon bars (SF24) are not Lanczos-upscaled; that step destroys the labels.
