"""
Radial efficiency charts for Generation streams-System {1,2,3}.xlsx.

Each figure: (a) energy efficiency (%), (b) exergy efficiency (%).
Style matches subsystem_energy_efficiency_radial.py.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

from config import DATA_DIR, chart_path

DATA_FILES = [
    DATA_DIR / "Generation streams-System 1.xlsx",
    DATA_DIR / "Generation streams-System 2.xlsx",
    DATA_DIR / "Generation streams-System 3.xlsx",
]

PALETTE = [
    "#e91e8c",
    "#7b2cbf",
    "#2563eb",
    "#0891b2",
    "#16a34a",
    "#84cc16",
    "#f59e0b",
    "#ea580c",
    "#dc2626",
    "#be123c",
]


def _load_system(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, header=0)
    df.columns = ["Stream", "Energy efficiency (%)", "Exergy efficiency (%)"]
    df["Stream"] = df["Stream"].astype(str).str.strip()
    for col in ("Energy efficiency (%)", "Exergy efficiency (%)"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["Stream"])


def _nice_axis_max(vmax: float) -> float:
    """Round axis ceiling upward to a clean percentage (≤ 100)."""
    if vmax <= 0:
        return 100.0
    for step in (5, 10, 20):
        ceiling = np.ceil(vmax / step) * step
        if ceiling >= vmax + 0.5 * step or ceiling >= 100:
            return float(min(100.0, max(ceiling, vmax)))
    return 100.0


def _short_name(name: str) -> str:
    """Use parenthetical abbreviation only (no brackets); else keep full name."""
    match = re.search(r"\(([^)]+)\)", name)
    return match.group(1).strip() if match else name.strip()


def _draw_radial(
    ax,
    ax_lab,
    names: list[str],
    values: list[float],
    *,
    fs_label: float,
    fs_tick: float,
    panel: str,
) -> None:
    order = np.argsort(values)
    names = [_short_name(names[i]) for i in order]
    values = [float(values[i]) for i in order]
    n = len(names)
    colors = (
        PALETTE[-n:]
        if n <= len(PALETTE)
        else [plt.cm.turbo(0.15 + 0.75 * i / max(n - 1, 1)) for i in range(n)]
    )

    vmax = max(values)
    axis_max = _nice_axis_max(vmax)

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    bar_height = 0.72
    gap = 0.28
    r0 = 1.2

    for i, (_name, val, color) in enumerate(zip(names, values, colors)):
        bottom = r0 + i * (bar_height + gap)
        width = val / axis_max * 2 * np.pi
        ax.bar(
            x=0.0,
            height=bar_height,
            width=width,
            bottom=bottom,
            color=to_rgba(color, 0.92),
            edgecolor="white",
            linewidth=1.0,
            align="edge",
            zorder=3,
        )
        ax.bar(
            x=0.0,
            height=bar_height,
            width=2 * np.pi,
            bottom=bottom,
            color=to_rgba("#d4d4d4", 0.35),
            edgecolor="none",
            align="edge",
            zorder=1,
        )

    r_outer = r0 + n * (bar_height + gap) - gap / 2

    for i in range(n):
        r_mid = r0 + i * (bar_height + gap) + bar_height / 2
        ax.plot(
            np.linspace(0, 2 * np.pi, 360),
            np.full(360, r_mid),
            color="#c0c0c0",
            lw=0.4,
            zorder=0,
        )

    # Percentage ticks on outer rim (every 5%; label every tick)
    tick_step = 5.0
    n_ticks = int(round(axis_max / tick_step))
    for k in range(n_ticks + 1):
        pct = k * tick_step
        if pct > axis_max + 1e-9:
            break
        # Full-circle value lands on 12:00 with "0" — skip duplicate max label/tick
        if pct >= axis_max - 1e-9:
            continue
        theta = pct / axis_max * 2 * np.pi
        major = k % 2 == 0
        tick_len = 0.18 if major else 0.10
        ax.plot(
            [theta, theta],
            [r_outer, r_outer + tick_len],
            color="black",
            lw=0.85 if major else 0.65,
            solid_capstyle="butt",
            zorder=5,
        )
        # 12:00 ("0") sits farther out; all other rim numbers keep the usual radius
        if pct == 0:
            r_label = r_outer + 1.05
        else:
            r_label = r_outer + (0.68 if major else 0.58)
        ax.text(
            theta,
            r_label,
            f"{pct:.0f}",
            ha="center",
            va="center",
            fontsize=fs_tick * (1.0 if major else 0.92),
            rotation=-np.degrees(theta),
            rotation_mode="anchor",
            zorder=6,
        )

    ax.plot(
        np.linspace(0, 2 * np.pi * (vmax / axis_max) * 1.02, 400),
        np.full(400, r_outer),
        color="black",
        lw=1.1,
        zorder=4,
    )

    ax.set_ylim(0, r_outer + 1.55)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)
    ax.set_title(f"({panel})", fontsize=fs_label * 1.2, fontweight="bold", pad=22)

    # Stream labels tucked just left of each bar start
    renderer = ax.figure.canvas.get_renderer()
    pad_to_bar = 0.012  # axes-fraction gap between label end and bar
    for i, (name, val) in enumerate(zip(names, values)):
        r_mid = r0 + i * (bar_height + gap) + bar_height / 2
        xd, yd = ax.transData.transform((0.0, r_mid))
        x_start, ya = ax.transAxes.inverted().transform((xd, yd))
        label = f"{name} - {int(round(val))}%"

        tmp = ax_lab.text(0.0, ya, label, ha="left", va="center", fontsize=fs_label)
        bb = tmp.get_window_extent(renderer=renderer)
        x0, _ = ax_lab.transData.inverted().transform((bb.x0, bb.y0))
        x1, _ = ax_lab.transData.inverted().transform((bb.x1, bb.y1))
        tmp.remove()
        text_width = abs(x1 - x0)
        new_left = x_start - text_width - pad_to_bar

        ax_lab.text(
            new_left,
            ya,
            label,
            ha="left",
            va="center",
            fontsize=fs_label,
            color="#1a1a1a",
            clip_on=False,
        )


def plot_system_figure(df: pd.DataFrame, output_path: str, dpi: int = 600) -> None:
    names = df["Stream"].tolist()
    energy = df["Energy efficiency (%)"].tolist()
    exergy = df["Exergy efficiency (%)"].tolist()

    # Side-by-side dials with a modest gap (no overlap)
    fig = plt.figure(figsize=(14.0, 8.4))
    fs_label = 15.0
    fs_tick = 14.5

    gs = fig.add_gridspec(
        1,
        2,
        left=0.02,
        right=0.99,
        top=0.93,
        bottom=0.03,
        wspace=0.10,
    )
    ax_a = fig.add_subplot(gs[0, 0], polar=True)
    ax_b = fig.add_subplot(gs[0, 1], polar=True)
    fig.canvas.draw()

    lab_a = fig.add_axes(ax_a.get_position(), frameon=False)
    lab_b = fig.add_axes(ax_b.get_position(), frameon=False)
    for lab in (lab_a, lab_b):
        lab.set_xlim(0, 1)
        lab.set_ylim(0, 1)
        lab.axis("off")

    _draw_radial(
        ax_a, lab_a, names, energy, fs_label=fs_label, fs_tick=fs_tick, panel="a"
    )
    _draw_radial(
        ax_b, lab_b, names, exergy, fs_label=fs_label, fs_tick=fs_tick, panel="b"
    )

    fig.savefig(
        output_path, dpi=dpi, bbox_inches="tight", facecolor="white", pad_inches=0.15
    )
    plt.close(fig)
    print(f"[OK] Saved {output_path}")


def main() -> None:
    for path in DATA_FILES:
        if not path.exists():
            raise FileNotFoundError(path)
        df = _load_system(path)
        stem = path.stem.replace(" ", "_").replace("-", "_").lower()
        out = chart_path(f"{stem}_efficiency_radial.png")
        plot_system_figure(df, out)


if __name__ == "__main__":
    main()
