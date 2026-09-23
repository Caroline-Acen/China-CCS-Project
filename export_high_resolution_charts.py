"""
Export manuscript-named copies into charts/high_resolution/.

PNGs are copied if both edges are already >= MIN_PX; otherwise they are
upscaled with Lanczos. Dense carbon-contribution grids are copied as-is.
Sankey HTML is copied without conversion.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from PIL import Image

from config import CHARTS_DIR, ensure_charts_dir

MIN_PX = 3000
HIGH_RES_DIRNAME = "high_resolution"

# dest in high_resolution/ -> source filename in charts/
PAPER_PNGS: dict[str, str] = {
    "figure1.png": "EOR+DSA pipeline avg.png",
    "figure2.png": "EOR+DSA tank avg.png",
    "figure3.png": "Cost difference (EOR+DSA tank avg - EOR+DSA pipeline avg).png",
    "figure4a.png": "seismic_risk_factor_map.png",
    "figure4b.png": "risk_storage_bubble.png",
    "figure5a.png": "storage_industry_map.png",
    "figure5b.png": "sink_suitability_map_panels.png",
    "figure6a.png": "storage_comparison_dsa.png",
    "figure6b.png": "storage_comparison_eor.png",
    "figure6c.png": "injection_unabated.png",
    "figure9a.png": "gasfield_sink_combined_gate_map.png",
    "figure9b.png": "hydrogen_hub_cost_supply_curve.png",
    "figure10a.png": "Carbon_DSA_bubble_chart.png",
    "figure10b.png": "Carbon_EOR_bubble_chart.png",
    "figure11.png": "figure11.png",
    "sup_figure1.png": "pga_multi_panel.png",
    "sup_figure2.png": "seismic_faults_by_name_map_v2.png",
    "sup_figure3.png": "vs30_map.png",
    "sup_figure5.png": "storage_gas_map.png",
    "sup_figure6.png": "DSA pipeline min.png",
    "sup_figure7.png": "EOR pipeline min.png",
    "sup_figure8.png": "DSA tank min.png",
    "sup_figure9.png": "EOR tank min.png",
    "sup_figure10.png": "DSA pipeline avg.png",
    "sup_figure11.png": "EOR pipeline avg.png",
    "sup_figure12.png": "DSA tank avg.png",
    "sup_figure13.png": "EOR tank avg.png",
    "sup_figure14.png": "DSA pipeline max.png",
    "sup_figure15.png": "EOR pipeline max.png",
    "sup_figure16.png": "DSA tank max.png",
    "sup_figure17.png": "EOR tank max.png",
    "sup_figure18.png": "injection_maps.png",
    "sup_figure19.png": "historical_seismicity_map.png",
    "sup_figure20.png": "sensitivity_pga_2.png",
    "sup_figure21.png": "sensitivity_amplification_1.png",
    "sup_figure22.png": "sensitivity_fault_lambda_2.png",
    "sup_figure23.png": "sensitivity_weight_abc.png",
}

PAPER_HTML: dict[str, str] = {
    "figure7.html": "DSA_Sankey_Plot.html",
    "figure8.html": "EOR_Sankey_Plot.html",
}

# Copy without Lanczos (labels degrade if upscaled).
COPY_AS_IS: dict[str, str] = {
    "sup_figure24a.png": "sup_figure24a.png",
    "sup_figure24b.png": "sup_figure24b.png",
}

OPTIONAL_PNGS: dict[str, str] = {
    "figure11m.png": "figure11m.png",
    "figure11_mean.png": "figure11_mean.png",
    "figure11_median.png": "figure11_median.png",
}

# Same map as main-text Figure 5a.
DUPLICATE_FROM_DEST: dict[str, str] = {
    "sup_figure4.png": "figure5a.png",
}

KEEP_EXTRA = {
    "figure11m.png",
    "figure11_mean.png",
    "figure11_median.png",
}

EXTRAS_TO_REMOVE = {
    "risk_storage_bubble_narrow.png",
}

Image.MAX_IMAGE_PIXELS = None


def high_res_dir() -> Path:
    out = ensure_charts_dir() / HIGH_RES_DIRNAME
    out.mkdir(parents=True, exist_ok=True)
    return out


def target_size(width: int, height: int, min_px: int = MIN_PX) -> tuple[int, int]:
    scale = max(min_px / width, min_px / height, 1.0)
    return (
        max(min_px, int(round(width * scale))),
        max(min_px, int(round(height * scale))),
    )


def export_one(src: Path, dest: Path, min_px: int = MIN_PX, *, as_is: bool = False) -> str:
    if as_is:
        shutil.copy2(src, dest)
        with Image.open(src) as im:
            return f"copy as-is {im.size[0]}x{im.size[1]}"

    with Image.open(src) as im:
        w, h = im.size
        tw, th = target_size(w, h, min_px)
        if (tw, th) == (w, h):
            shutil.copy2(src, dest)
            return f"copy {w}x{h}"

        if im.mode not in ("RGB", "RGBA"):
            im = im.convert("RGBA" if "A" in im.getbands() else "RGB")
        out = im.resize((tw, th), Image.Resampling.LANCZOS)
        out.save(dest, format="PNG", optimize=True)
        return f"upscale {w}x{h} -> {tw}x{th}"


def _export_mapped(mapping: dict[str, str], *, as_is: bool = False) -> list[str]:
    charts = ensure_charts_dir()
    out_dir = high_res_dir()
    written: list[str] = []
    for dest_name, src_name in mapping.items():
        src = charts / src_name
        dest = out_dir / dest_name
        if not src.exists():
            if dest.exists():
                print(f"  KEEP existing {dest_name} (missing source {src_name})")
            else:
                print(f"  SKIP {dest_name}: missing {src_name}")
            continue
        action = export_one(src, dest, as_is=as_is)
        print(f"  {dest_name}: {action}  [{src_name}]")
        written.append(dest_name)
    return written


def main() -> None:
    out_dir = high_res_dir()
    print(f"Exporting manuscript figures to {out_dir} (min {MIN_PX}x{MIN_PX})")

    _export_mapped(PAPER_PNGS)
    _export_mapped(COPY_AS_IS, as_is=True)
    _export_mapped(OPTIONAL_PNGS)

    charts = ensure_charts_dir()
    for dest_name, src_name in PAPER_HTML.items():
        src = charts / src_name
        dest = out_dir / dest_name
        if not src.exists():
            print(f"  SKIP {dest_name}: missing {src_name}")
            continue
        shutil.copy2(src, dest)
        print(f"  {dest_name}: copy HTML  [{src_name}]")

    for dest_name, from_dest in DUPLICATE_FROM_DEST.items():
        src = out_dir / from_dest
        dest = out_dir / dest_name
        if not src.exists():
            print(f"  SKIP {dest_name}: missing {from_dest}")
            continue
        shutil.copy2(src, dest)
        print(f"  {dest_name}: copy of {from_dest}")

    keep = (
        set(PAPER_PNGS)
        | set(PAPER_HTML)
        | set(COPY_AS_IS)
        | set(DUPLICATE_FROM_DEST)
        | set(OPTIONAL_PNGS)
        | KEEP_EXTRA
    )
    for leftover in list(out_dir.iterdir()):
        if leftover.is_file() and leftover.name not in keep:
            leftover.unlink()
            print(f"  removed extra {leftover.name}")

    for extra in EXTRAS_TO_REMOVE:
        stale = out_dir / extra
        if stale.exists():
            stale.unlink()
            print(f"  removed extra {extra}")

    written = sorted(p.name for p in out_dir.iterdir() if p.is_file())
    print(f"Done. {len(written)} files in {out_dir}")


if __name__ == "__main__":
    main()
