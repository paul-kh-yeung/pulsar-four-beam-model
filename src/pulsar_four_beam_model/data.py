"""Input/output helpers for binned pulsar phaseograms."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
from scipy.special import gammaln

Array = np.ndarray


@dataclass(frozen=True)
class Phaseogram:
    """Binned phaseogram with optional per-bin uncertainties.

    The likelihood implemented in this package follows the weighted Poisson
    form used in the original research script.
    """

    phase: Array
    counts: Array
    count_errors: Array
    weights: Array
    log_factorials: Array
    source_path: str | None = None

    @property
    def n_bins(self) -> int:
        return int(self.phase.size)


def load_phaseogram(
    path: str | Path,
    phase_window: Tuple[float, float] = (0.80001, 1.80001),
    phase_shift: float = -1.0,
) -> Phaseogram:
    """Load a text phaseogram file.

    The expected columns are ``phase count [count_error]``. If the third column
    is absent, Poisson errors ``sqrt(count)`` are used. The default phase window
    and shift reproduce the original Vela scripts, where the interval
    0.80001--1.80001 is selected and shifted by -1.
    """
    path = Path(path)
    data = np.loadtxt(path, dtype=float)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"{path} must contain at least two columns: phase and counts")

    phase_raw = data[:, 0]
    mask = (phase_raw > phase_window[0]) & (phase_raw < phase_window[1])
    data = data[mask]
    if data.size == 0:
        raise ValueError(f"No bins selected from {path}; check phase_window={phase_window}")

    phase = data[:, 0] + phase_shift
    order = np.argsort(phase)
    phase = phase[order]
    counts = data[:, 1][order]
    if data.shape[1] >= 3:
        count_errors = data[:, 2][order]
    else:
        count_errors = np.sqrt(np.clip(counts, 0.0, np.inf))

    weights = np.round(counts / np.clip(count_errors, 1e-300, np.inf) ** 2)
    weights[~(weights > 0)] = 1.0
    log_factorials = gammaln(counts * weights + 1.0)
    return Phaseogram(
        phase=phase,
        counts=counts,
        count_errors=count_errors,
        weights=weights,
        log_factorials=log_factorials,
        source_path=str(path),
    )


def load_many(paths: Iterable[str | Path]) -> dict[str, Phaseogram]:
    """Load several phaseogram files and return them by stem."""
    result: dict[str, Phaseogram] = {}
    for path in paths:
        path = Path(path)
        result[path.stem] = load_phaseogram(path)
    return result


def estimate_off_pulse_systematic(counts: Array, phase: Array, off_phase: Tuple[float, float] = (-0.2, 0.08)) -> float:
    """Estimate the fractional off-pulse systematic used for plotting checks.

    This is included for reproducibility of the original diagnostic plots. It is
    not used by the default likelihood.
    """
    mask = (phase >= off_phase[0]) & (phase <= off_phase[1])
    off_counts = np.asarray(counts, dtype=float)[mask]
    if off_counts.size == 0:
        raise ValueError("No off-pulse bins found; check off_phase")
    variance = float(np.std(off_counts) ** 2)
    mean = float(np.mean(off_counts))
    if mean <= 0:
        return 0.0
    excess_variance = variance - min(variance, mean)
    return float(np.sqrt(max(excess_variance, 0.0)) / mean)
