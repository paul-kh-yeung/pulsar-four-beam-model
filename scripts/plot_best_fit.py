#!/usr/bin/env python
"""Plot a saved four-beam solution.

Example
-------
python scripts/plot_best_fit.py configs/vela_v1_97.json --output figures/vela_best_fit.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pulsar_four_beam_model.config import read_json
from pulsar_four_beam_model.data import load_phaseogram
from pulsar_four_beam_model.likelihood import joint_deviance
from pulsar_four_beam_model.model import band_model_from_flat_parameters, geometry_from_flat_parameters
from pulsar_four_beam_model.plotting import plot_all_bands
from pulsar_four_beam_model.derived import (
    compute_derived_geometry,
    format_derived_geometry,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument("--data-root", default=".", help="Root directory prepended to data paths in the config")
    parser.add_argument("--output", default="figures/four_beam_best_fit.png", help="Output figure path")
    args = parser.parse_args()

    config_path = Path(args.config)
    cfg = read_json(config_path)
    params = cfg["parameters"]
    geometry = geometry_from_flat_parameters(params)
    band_order = cfg.get("band_order", list(cfg["bands"].keys()))
    data_root = Path(args.data_root)

    phaseograms = {}
    band_models = {}
    titles = {}
    for band in band_order:
        band_cfg = cfg["bands"][band]
        phaseograms[band] = load_phaseogram(data_root / band_cfg["file"])
        band_models[band] = band_model_from_flat_parameters(params, band)
        titles[band] = band_cfg.get("label", band)

    deviance = joint_deviance(phaseograms, geometry, band_models)
    print(f"Joint -2 log L = {deviance:.6g}")
    print("Axis phases:")
    for name, value in geometry.axis_phases().items():
        print(f"  {name}: {value:.8f}")

    derived = compute_derived_geometry(params)
    print("\nDerived geometric quantities")
    print("----------------------------")
    print(format_derived_geometry(derived))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    plot_all_bands(phaseograms, geometry, band_models, band_order, titles=titles, output=output)
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
