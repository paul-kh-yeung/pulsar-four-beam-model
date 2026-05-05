"""Reference implementation of the pulsar four-beam model."""

from .model import (
    BeamPairGeometry,
    BeamShape,
    EnergyBandModel,
    FourBeamGeometry,
    band_model_from_flat_parameters,
    geometry_from_flat_parameters,
    simulate_four_beam_band,
)

from pulsar_four_beam_model.derived import (
    DerivedGeometry,
    compute_derived_geometry,
    format_derived_geometry,
)

__all__ = [
    "BeamPairGeometry",
    "BeamShape",
    "EnergyBandModel",
    "FourBeamGeometry",
    "band_model_from_flat_parameters",
    "geometry_from_flat_parameters",
    "simulate_four_beam_band",
    "DerivedGeometry",
    "compute_derived_geometry",
    "format_derived_geometry",
]
