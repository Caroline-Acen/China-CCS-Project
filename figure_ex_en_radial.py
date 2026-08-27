"""
Radial bar charts for subsystem energy / exergy efficiency
(Ex-En-Subsystem.xlsx). Values shown as fractions, not percentages.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

from config import DATA_DIR, chart_path, manuscript_font_bundle

DATA_XLSX = DATA_DIR / "Ex-En-Subsystem.xlsx"

# Rainbow-ish palette (inner → outer), similar to the reference chart
PALETTE = [
    "#e91e8c",  # magenta
    "#7b2cbf",  # purple
    "#2563eb",  # blue
    "#0891b2",  # cyan
    "#16a34a",  # green
    "#84cc16",  # lime
    "#f59e0b",  # amber
    "#ea580c",  # orange
    "#dc2626",  # red
]


def _load_efficiencies() -> pd.DataFrame:
    df = pd.read_excel(DATA_XLSX)
    df.columns = [str(c).strip() for c in df.columns]
    if "Subsystem" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Subsystem"})
    df["Subsystem"] = df["Subsystem"].astype(str).str.strip()
    return df


def plot_radial_efficiency(
    names: list[str],
    values: list[float],
    output_path: str,
    dpi: int = 600,
) -> None:
    """Concentric radial bars; outermost = largest value (sorted ascending)."""
    order = np.argsort(values)
    names = [names[i] for i in order]
    values = [float(values[i]) for i in order]
    n = len(names)
    colors = PALETTE[-n:] if n <= len(PALETTE) else [
        plt.cm.turbo(0.15 + 0.75 * i / max(n - 1, 1)) for i in range(n)
    ]

    vmax = max(values)
    axis_max = min(1.0, np.ceil(vmax * 20) / 20 + 0.05)
    axis_max = max(axis_max, vmax + 0.02)

    fig = plt.figure(figsize=(8.5, 8.5))
    mf = manuscript_font_bundle(fig.get_figwidth())
    fs_label = mf["label"] * 1.25
    fs_tick = mf["tick"] * 1.4

    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)  # clockwise

    bar_height = 0.72
    gap = 0.28
    r0 = 1.2  # hollow center

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

    # Outer fraction axis
    tick_step = 0.05
    n_ticks = int(round(axis_max / tick_step))
    for k in range(n_ticks + 1):
        frac = k * tick_step
        if frac > axis_max + 1e-9:
            break
        theta = frac / axis_max * 2 * np.pi
        tick_len = 0.12 if k % 2 == 0 else 0.06
        ax.plot(
            [theta, theta],
            [r_outer, r_outer + tick_len],
            color="black",
            lw=0.7,
            solid_capstyle="butt",
            zorder=5,
        )
        if k % 2 == 0:
            label = "0" if frac == 0 else f"{frac:.2f}".rstrip("0").rstrip(".")
            ax.text(
                theta,
                r_outer + 0.42,
                label,
                ha="center",
                va="center",
                fontsize=fs_tick,
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

    ax.set_ylim(0, r_outer + 0.85)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)

    # Labels at bar starts (left of 12 o'clock), no parentheses
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_lab = fig.add_axes(ax.get_position(), frameon=False)
    ax_lab.set_xlim(0, 1)
    ax_lab.set_ylim(0, 1)
    ax_lab.axis("off")

    # Previous left anchor; move each label so the gap to the bar start is halved
    prev_label_x = 0.28 + 10.0 / ax.get_window_extent().width

    for i, (name, val) in enumerate(zip(names, values)):
        r_mid = r0 + i * (bar_height + gap) + bar_height / 2
        xd, yd = ax.transData.transform((0.0, r_mid))
        x_start, ya = ax.transAxes.inverted().transform((xd, yd))
        label = f"{name} - {val:.2f}"

        # Measure text width in axes coordinates
        tmp = ax_lab.text(
            prev_label_x, ya, label, ha="left", va="center", fontsize=fs_label
        )
        bb = tmp.get_window_extent(renderer=renderer)
        x0, _ = ax_lab.transData.inverted().transform((bb.x0, bb.y0))
        x1, _ = ax_lab.transData.inverted().transform((bb.x1, bb.y1))
        tmp.remove()
        text_right = max(x0, x1)
        gap_to_bar = x_start - text_right
        new_left = prev_label_x + 0.75 * gap_to_bar

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

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white", pad_inches=0.15)
    plt.close(fig)
    print(f"[OK] Saved {output_path}")


def main() -> None:
    df = _load_efficiencies()

    energy = df[df["Subsystem"] != "AGMD"]
    plot_radial_efficiency(
        energy["Subsystem"].tolist(),
        energy["Energy efficiency"].tolist(),
        output_path=chart_path("subsystem_energy_efficiency_radial.png"),
    )
    plot_radial_efficiency(
        df["Subsystem"].tolist(),
        df["Exergy efficiency"].tolist(),
        output_path=chart_path("subsystem_exergy_efficiency_radial.png"),
    )


if __name__ == "__main__":
    main()
