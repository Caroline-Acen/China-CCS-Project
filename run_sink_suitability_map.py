from pathlib import Path

from figure_sink_suitability_distribution import (
    chart_path,
    ensure_charts_dir,
    plot_map_panels,
)


def main() -> None:
    ensure_charts_dir()
    output_path = chart_path("sink_suitability_map_panels.png")
    plot_map_panels(
        [2025, 2030, 2035, 2040, 2045, 2050],
        output_path,
        combined_path=None,
        categories=None,
    )


if __name__ == "__main__":
    main()
