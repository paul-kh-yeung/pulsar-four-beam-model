#!/usr/bin/env python
"""Plot beam orientation psi and inverse Doppler factor epsilon.

Example
-------
python scripts/plot_psi_epsilon.py configs/vela_v1_97.json \
  --output figures/vela_psi_epsilon.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pulsar_four_beam_model.config import read_json
from pulsar_four_beam_model.data import load_phaseogram
from pulsar_four_beam_model.model import (
    band_model_from_flat_parameters,
    geometry_from_flat_parameters,
)
from pulsar_four_beam_model.plotting import plot_psi_epsilon


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument(
        "--data-root",
        default=".",
        help="Root directory prepended to data paths in the config",
    )
    parser.add_argument(
        "--reference-band",
        default=None,
        help=(
            "Energy-band prefix used only to provide the phase grid. "
            "Default: first band in band_order."
        ),
    )
    parser.add_argument(
        "--output",
        default="figures/psi_epsilon.png",
        help="Output figure path",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Figure title. Default: config title or empty string.",
    )
    args = parser.parse_args()

    cfg = read_json(args.config)
    params = cfg["parameters"]

    geometry = geometry_from_flat_parameters(params)
    band_order = cfg.get("band_order", list(cfg["bands"].keys()))
    reference_band = args.reference_band or band_order[0]

    if reference_band not in cfg["bands"]:
        raise KeyError(
            f"Reference band {reference_band!r} not found in config bands: "
            f"{list(cfg['bands'].keys())}"
        )

    data_root = Path(args.data_root)
    band_cfg = cfg["bands"][reference_band]
    phaseogram = load_phaseogram(data_root / band_cfg["file"])
    band_model = band_model_from_flat_parameters(params, reference_band)

    title = args.title
    if title is None:
        title = cfg.get("name", "")

    output = Path(args.output)
    plot_psi_epsilon(
        phaseogram=phaseogram,
        geometry=geometry,
        band_model=band_model,
        title=title,
        output=output,
    )

    print(f"Reference band: {reference_band}")
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
