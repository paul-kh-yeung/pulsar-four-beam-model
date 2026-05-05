#!/usr/bin/env python
"""Evaluate the joint likelihood of a saved four-beam solution."""

from __future__ import annotations

import argparse
from pathlib import Path

from pulsar_four_beam_model.config import read_json
from pulsar_four_beam_model.data import load_phaseogram
from pulsar_four_beam_model.likelihood import joint_deviance
from pulsar_four_beam_model.model import band_model_from_flat_parameters, geometry_from_flat_parameters
from pulsar_four_beam_model.derived import (
    compute_derived_geometry,
    format_derived_geometry,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument("--data-root", default=".", help="Root directory prepended to data paths in the config")
    args = parser.parse_args()

    cfg = read_json(args.config)
    derived = compute_derived_geometry(cfg["parameters"])

    print("\nDerived geometric quantities")
    print("----------------------------")
    print(format_derived_geometry(derived))
    params = cfg["parameters"]
    geometry = geometry_from_flat_parameters(params)
    phaseograms = {}
    band_models = {}
    for band in cfg.get("band_order", list(cfg["bands"].keys())):
        phaseograms[band] = load_phaseogram(Path(args.data_root) / cfg["bands"][band]["file"])
        band_models[band] = band_model_from_flat_parameters(params, band)
    print(joint_deviance(phaseograms, geometry, band_models))


if __name__ == "__main__":
    main()
