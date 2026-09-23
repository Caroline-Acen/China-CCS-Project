"""Figure 11 variant with a 250 km multiplier on all transport values."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import figure11
from config import chart_path
from export_high_resolution_charts import export_one, high_res_dir

TRANSPORT_DISTANCE_KM = 250.0


def main() -> None:
    saved = deepcopy(figure11.DATA)
    try:
        # Scale every country's transport. China in figure11.DATA is already
        # USD/ton (×250 km); leave it. USA and EU are still per-ton totals.
        for country, values in figure11.DATA.items():
            if country == "China":
                continue
            lo, hi = values["Transportation"]
            values["Transportation"] = (
                lo * TRANSPORT_DISTANCE_KM,
                hi * TRANSPORT_DISTANCE_KM,
            )

        output = Path(chart_path("figure11m.png"))
        figure11.plot_figure11(str(output))
        high_res_output = high_res_dir() / output.name
        print(export_one(output, high_res_output))
    finally:
        figure11.DATA.clear()
        figure11.DATA.update(saved)


if __name__ == "__main__":
    main()
