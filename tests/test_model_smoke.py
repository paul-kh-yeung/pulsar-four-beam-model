import numpy as np

from pulsar_four_beam_model.model import (
    BeamPairGeometry,
    BeamShape,
    EnergyBandModel,
    FourBeamGeometry,
    simulate_four_beam_band,
)


def test_four_beam_simulation_smoke():
    phase = np.linspace(-0.2, 0.8, 64)
    geometry = FourBeamGeometry(
        theta_A_deg=54.85,
        lower=BeamPairGeometry(
            r_sin_colat_over_rlc=1.0,
            phase_separation=0.43,
            t0_phase=0.29,
            gamma_lor=10.0,
            bulk_zenith_deg=146.0,
            bulk_azimuth_deg=83.0,
            beam_zenith_deg=123.0,
            beam_azimuth_deg=88.0,
        ),
        higher=BeamPairGeometry(
            r_sin_colat_over_rlc=0.6,
            phase_separation=0.25,
            t0_phase=0.39,
            gamma_lor=1.4,
            bulk_zenith_deg=31.0,
            bulk_azimuth_deg=37.0,
            beam_zenith_deg=34.0,
            beam_azimuth_deg=37.0,
        ),
    )
    band = EnergyBandModel(
        lower_shape=BeamShape(log10_A=1.0, eta=1.0, psi_c_deg=20.0, beta=2.0),
        higher_shape=BeamShape(log10_A=1.0, eta=1.0, psi_c_deg=10.0, beta=2.0),
        background=5.0,
    )
    sim = simulate_four_beam_band(phase, geometry, band)
    assert sim.model.shape == phase.shape
    assert np.all(np.isfinite(sim.model))
    assert np.all(sim.model > 0)
