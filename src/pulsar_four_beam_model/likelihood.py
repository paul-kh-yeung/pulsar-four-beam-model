"""Likelihood utilities for four-beam phaseogram fitting."""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .data import Phaseogram
from .model import EnergyBandModel, FourBeamGeometry, simulate_four_beam_band


def weighted_poisson_deviance(phaseogram: Phaseogram, model_counts: np.ndarray) -> float:
    """Return ``-2 log L`` for the weighted Poisson likelihood.

    This matches the form used by the original fitting script:
    ``lnP = -lambda*w + log(lambda*w)*count*w - log((count*w)!)``.
    """
    model = np.asarray(model_counts, dtype=float)
    model = np.clip(model, 1e-300, np.inf)
    w = phaseogram.weights
    c = phaseogram.counts
    logp = -model * w + np.log(model * w) * c * w - phaseogram.log_factorials
    return float(-2.0 * np.sum(logp))


def joint_deviance(
    phaseograms: Mapping[str, Phaseogram],
    geometry: FourBeamGeometry,
    band_models: Mapping[str, EnergyBandModel],
    interpolation_grid_size: int = 2200,
) -> float:
    """Return the joint ``-2 log L`` over all supplied energy bands."""
    total = 0.0
    for band_name, phaseogram in phaseograms.items():
        if band_name not in band_models:
            raise KeyError(f"Missing model parameters for energy band {band_name!r}")
        simulation = simulate_four_beam_band(
            phaseogram.phase,
            geometry,
            band_models[band_name],
            interpolation_grid_size=interpolation_grid_size,
        )
        total += weighted_poisson_deviance(phaseogram, simulation.model)
    return float(total)
