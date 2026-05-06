"""Minuit fitting helpers for the four-beam model.

The public API fits the 53-parameter compact model used in the accepted paper:
17 geometry/kinematic parameters plus 9 band-dependent parameters per energy
band (for four bands: 17 + 4*9 = 53).
"""

from __future__ import annotations

from typing import Callable, Mapping, Sequence

from iminuit import Minuit

from .data import Phaseogram
from .likelihood import joint_deviance
from .model import band_model_from_flat_parameters, geometry_from_flat_parameters

GEOMETRY_PARAMETER_NAMES = [
    "angleA_deg",
    "d1sin",
    "P2minusP1",
    "t___0",
    "lor",
    "thetaM_deg",
    "thetaN_deg",
    "thetaC_deg",
    "thetaQ_deg",
    "D3SIN",
    "P4minusP3",
    "T___O",
    "LOR",
    "phiM_deg",
    "phiN_deg",
    "phiC_deg",
    "phiQ_deg",
]

BAND_PARAMETER_SUFFIXES = [
    "LogNorm0",
    "Eta1",
    "Psi_c0",
    "Beta",
    "LogAmpO",
    "eTA1",
    "pSI_wO",
    "bETA",
    "Const",
]


def active_parameter_names(bands: Sequence[str] = ("BM", "LO", "HI", "TP")) -> list[str]:
    """Return the 53 active parameter names used by the compact benchmark model."""
    names = list(GEOMETRY_PARAMETER_NAMES)
    for band in bands:
        names.extend(f"{band}_{suffix}" for suffix in BAND_PARAMETER_SUFFIXES)
    return names


def make_joint_objective(
    phaseograms: Mapping[str, Phaseogram],
    base_parameters: Mapping[str, float],
    bands: Sequence[str] = ("BM", "LO", "HI", "TP"),
    interpolation_grid_size: int = 2200,
) -> Callable[..., float]:
    """Create a Minuit-compatible objective function.

    A small dynamic wrapper is used so that Minuit sees a normal Python
    signature with named scalar parameters.
    """
    parameter_names = active_parameter_names(bands)
    missing = [name for name in parameter_names if name not in base_parameters]
    if missing:
        raise KeyError(f"Missing starting values for active parameters: {missing}")

    def evaluate(values: Mapping[str, float]) -> float:
        params = dict(base_parameters)
        params.update(values)
        geometry = geometry_from_flat_parameters(params)
        band_models = {band: band_model_from_flat_parameters(params, band) for band in bands}
        return joint_deviance(
            {band: phaseograms[band] for band in bands},
            geometry,
            band_models,
            interpolation_grid_size=interpolation_grid_size,
        )

    args = ", ".join(parameter_names)
    mapping = ", ".join(f"'{name}': {name}" for name in parameter_names)
    namespace = {"evaluate": evaluate}
    source = f"def objective({args}):\n    return evaluate({{{mapping}}})\n"
    exec(source, namespace)
    return namespace["objective"]


def build_minuit(
    phaseograms: Mapping[str, Phaseogram],
    starting_parameters: Mapping[str, float],
    bands: Sequence[str] = ("BM", "LO", "HI", "TP"),
    interpolation_grid_size: int = 2200,
) -> Minuit:
    """Build a Minuit object from phaseograms and starting parameters."""
    objective = make_joint_objective(
        phaseograms,
        starting_parameters,
        bands=bands,
        interpolation_grid_size=interpolation_grid_size,
    )
    start = {name: float(starting_parameters[name]) for name in active_parameter_names(bands)}
    minuit = Minuit(objective, **start)
    minuit.errordef = 1.0  # Objective is -2 log L.
    apply_default_limits(minuit, bands=bands)
    return minuit


def apply_default_limits(minuit: Minuit, bands: Sequence[str]) -> None:
    """Apply broad, conservative limits similar to the legacy fitting script."""
    minuit.limits["lor"] = (1.01, 91.0)
    minuit.limits["LOR"] = (1.01, 91.0)
    minuit.limits["d1sin"] = (0.3, 1.2)
    minuit.limits["D3SIN"] = (0.3, 1.2)
    minuit.limits["t___0"] = (-0.3, 0.7)
    minuit.limits["T___O"] = (0.0, 1.0)
    minuit.limits["P2minusP1"] = (0.2, 0.45)
    minuit.limits["P4minusP3"] = (0.2, 0.45)
    minuit.limits["thetaC_deg"] = (10.0, 170.0)
    minuit.limits["thetaQ_deg"] = (60.0, 335.0)
    minuit.limits["phiC_deg"] = (0.2, 170.0)
    minuit.limits["phiQ_deg"] = (-270.0, 90.0)
    minuit.limits["thetaM_deg"] = (90.0, 180.0)
    minuit.limits["thetaN_deg"] = (-90.0, 270.0)
    minuit.limits["phiM_deg"] = (0.0, 90.0)
    minuit.limits["phiN_deg"] = (-270.0, 90.0)
    minuit.limits["angleA_deg"] = (0.0, 90.0)
    for band in bands:
        minuit.limits[f"{band}_Beta"] = (0.01, 9.0)
        minuit.limits[f"{band}_bETA"] = (0.2, 9.0)
        minuit.limits[f"{band}_Psi_c0"] = (0.1, 180.0)
        minuit.limits[f"{band}_pSI_wO"] = (0.1, 180.0)
        minuit.limits[f"{band}_Const"] = (0.0, None)
