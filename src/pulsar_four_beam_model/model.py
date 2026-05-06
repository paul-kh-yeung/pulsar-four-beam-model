r"""Geometry-first four-beam model for gamma-ray pulsar phaseograms.

The notation follows Yeung & Saito (ApJ, accepted):

* ``theta_A``: inclination angle between the spin axis and the line of sight.
* ``d_sin_theta_B``: cylindrical radius of the lower-altitude site in units of
  the light-cylinder radius, :math:`d\sin\theta_B/R_{LC}`.
* ``D_sin_phi_B``: cylindrical radius of the higher-altitude site in units of
  the light-cylinder radius, :math:`D\sin\phi_B/R_{LC}`.
* ``gamma_lor``: effective bulk Lorentz factor.
* ``theta_M, theta_N`` or ``phi_M, phi_N``: zenith and azimuthal angles of the
  bulk-flow direction.
* ``theta_C, theta_Q`` or ``phi_C, phi_Q``: zenith and azimuthal angles of the
  beam axis.
* ``epsilon``: inverse Doppler factor, :math:`E_{bulk}/E_{det}`.
* ``psi``: beam-axis--line-of-sight angle.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, MutableMapping

import numpy as np

Array = np.ndarray
DEG_TO_RAD = np.pi / 180.0
RAD_TO_DEG = 180.0 / np.pi


def _as_float_array(values: Array | float) -> Array:
    """Return ``values`` as a NumPy floating array."""
    return np.asarray(values, dtype=float)


def _clip_for_arccos(values: Array) -> Array:
    """Avoid nan values caused by tiny floating-point excursions beyond [-1, 1]."""
    return np.clip(values, -1.0, 1.0)


def lorentz_beta(gamma_lor: float) -> float:
    """Return :math:`v/c` for a bulk Lorentz factor.

    Parameters
    ----------
    gamma_lor:
        Bulk Lorentz factor. Must be larger than 1.
    """
    if gamma_lor <= 1.0:
        raise ValueError("gamma_lor must be larger than 1")
    return float(np.sqrt(1.0 - 1.0 / gamma_lor**2))


@dataclass(frozen=True)
class BeamShape:
    """Band-dependent generalized-normal-distribution beam shape.

    This is the compact form used in the ApJ paper,

    ``C_i(psi, epsilon) = A_i epsilon**eta_i exp[-(psi/Psi_c_i)**beta_i]``.

    Extra legacy parameters are deliberately omitted here because the accepted
    four-beam implementation fixed those higher-order terms to zero.
    """

    log10_A: float
    eta: float
    psi_c_deg: float
    beta: float

    def count_rate(self, psi_deg: Array, epsilon: Array) -> Array:
        """Evaluate the beam contribution to a phaseogram."""
        psi = _as_float_array(psi_deg)
        eps = _as_float_array(epsilon)
        eps = np.clip(eps, 1e-300, np.inf)
        psi_c = max(float(self.psi_c_deg), 1e-12)
        exponent = -np.power(np.clip(psi / psi_c, 0.0, np.inf), self.beta)
        log_rate = np.log(10.0) * self.log10_A + self.eta * np.log(eps) + exponent
        return np.exp(np.clip(log_rate, -745.0, 709.0))


@dataclass(frozen=True)
class BeamPairGeometry:
    """Shared geometry and bulk motion for one north--south beam pair."""

    r_sin_colat_over_rlc: float
    phase_separation: float
    t0_phase: float
    gamma_lor: float
    bulk_zenith_deg: float
    bulk_azimuth_deg: float
    beam_zenith_deg: float
    beam_azimuth_deg: float

    @property
    def beta_bulk(self) -> float:
        """Return bulk speed in units of ``c``."""
        return lorentz_beta(self.gamma_lor)

    def colatitude_rad(self, theta_A_deg: float) -> float:
        """Infer the site colatitude from the fitted phase separation.

        This preserves the algebra used in the original fitting script, which
        follows the time-delay relation between the two antipodal beams.
        """
        theta_A = theta_A_deg * DEG_TO_RAD
        phase_delta = self.phase_separation - np.floor(self.phase_separation) - 0.5
        if np.isclose(phase_delta, 0.0):
            raise ValueError("phase_separation is too close to 0.5; colatitude is singular")
        argument = -self.r_sin_colat_over_rlc * np.cos(theta_A) / (np.pi * phase_delta)
        colat = float(np.arctan(argument))
        if colat < 0.0:
            colat += np.pi
        return colat

    def radius_over_rlc(self, theta_A_deg: float) -> float:
        """Return representative site radius in units of ``R_LC``."""
        colat = self.colatitude_rad(theta_A_deg)
        return float(self.r_sin_colat_over_rlc / np.sin(colat))

    def vertical_over_rlc(self, theta_A_deg: float) -> float:
        """Return representative site height in units of ``R_LC``."""
        colat = self.colatitude_rad(theta_A_deg)
        return float(self.radius_over_rlc(theta_A_deg) * np.cos(colat))

    def first_axis_phase(self, theta_A_deg: float) -> float:
        """Observed phase when the first beam axis is closest to the line of sight."""
        theta_A = theta_A_deg * DEG_TO_RAD
        q = self.beam_azimuth_deg * DEG_TO_RAD
        return float(
            self.t0_phase
            - self.beam_azimuth_deg / 360.0
            - self.r_sin_colat_over_rlc / (2.0 * np.pi) * (np.cos(-q) - 1.0) * np.sin(theta_A)
        )

    def second_axis_phase(self, theta_A_deg: float) -> float:
        """Observed phase when the antipodal beam axis is closest to the line of sight."""
        return self.first_axis_phase(theta_A_deg) + self.phase_separation


@dataclass(frozen=True)
class FourBeamGeometry:
    """Common geometry for the two lower- and higher-altitude beam pairs."""

    theta_A_deg: float
    lower: BeamPairGeometry
    higher: BeamPairGeometry

    def axis_phases(self) -> Dict[str, float]:
        """Return the four phases at which beam axes are closest to the line of sight."""
        return {
            "beam1": self.lower.first_axis_phase(self.theta_A_deg),
            "beam2": self.lower.second_axis_phase(self.theta_A_deg),
            "beam3": self.higher.first_axis_phase(self.theta_A_deg),
            "beam4": self.higher.second_axis_phase(self.theta_A_deg),
        }


@dataclass(frozen=True)
class EnergyBandModel:
    """Band-dependent model parameters for one phaseogram."""

    lower_shape: BeamShape
    higher_shape: BeamShape
    background: float


@dataclass(frozen=True)
class BeamPairSimulation:
    """Phase-dependent quantities for one beam pair."""

    colatitude_rad: float
    radius_over_rlc: float
    velocity_over_c: float
    psi_first_deg: Array
    psi_second_deg: Array
    epsilon_first: Array
    epsilon_second: Array
    counts_first: Array
    counts_second: Array


@dataclass(frozen=True)
class FourBeamSimulation:
    """Four-beam model prediction for one energy band."""

    phase: Array
    background: float
    beam1: Array
    beam2: Array
    beam3: Array
    beam4: Array
    model: Array
    lower_pair: BeamPairSimulation
    higher_pair: BeamPairSimulation

    @property
    def components(self) -> Mapping[str, Array]:
        """Return components as a name-to-array mapping."""
        return {
            "background": np.full_like(self.phase, self.background, dtype=float),
            "beam1": self.beam1,
            "beam2": self.beam2,
            "beam3": self.beam3,
            "beam4": self.beam4,
            "model": self.model,
        }


def simulate_beam_pair(
    phase: Array,
    theta_A_deg: float,
    pair: BeamPairGeometry,
    shape: BeamShape,
    interpolation_grid_size: int = 2200,
) -> BeamPairSimulation:
    """Simulate one north--south beam pair.

    The implementation intentionally mirrors the original research script, but
    exposes the variables with paper-style notation.
    """
    phase = _as_float_array(phase)
    theta_A = theta_A_deg * DEG_TO_RAD
    theta_C = pair.beam_zenith_deg * DEG_TO_RAD
    theta_Q = pair.beam_azimuth_deg * DEG_TO_RAD
    theta_M = pair.bulk_zenith_deg * DEG_TO_RAD
    theta_N = pair.bulk_azimuth_deg * DEG_TO_RAD

    colatitude = pair.colatitude_rad(theta_A_deg)
    radius = pair.radius_over_rlc(theta_A_deg)
    beta_bulk = pair.beta_bulk

    def psi_second(z: Array) -> Array:
        arg = -np.cos(z * 2.0 * np.pi + theta_Q) * np.sin(theta_C) * np.sin(theta_A) + np.cos(theta_C) * np.cos(theta_A)
        return RAD_TO_DEG * np.arccos(_clip_for_arccos(arg))

    def psi_first(z: Array) -> Array:
        arg = np.cos(z * 2.0 * np.pi + theta_Q) * np.sin(theta_C) * np.sin(theta_A) - np.cos(theta_C) * np.cos(theta_A)
        return RAD_TO_DEG * np.arccos(_clip_for_arccos(arg))

    first_axis_phase = pair.first_axis_phase(theta_A_deg)
    y0_det = -theta_Q / (2.0 * np.pi) - radius / (2.0 * np.pi) * (
        np.cos(-theta_Q) * np.sin(colatitude) * np.sin(theta_A) - np.cos(colatitude) * np.cos(theta_A)
    )

    grid = np.linspace(-1.1, 1.1, interpolation_grid_size)
    x_psr = grid
    y_psr = grid
    x_det = x_psr - radius / (2.0 * np.pi) * (
        -np.cos(x_psr * 2.0 * np.pi) * np.sin(colatitude) * np.sin(theta_A) + np.cos(colatitude) * np.cos(theta_A)
    )
    y_det = y_psr - radius / (2.0 * np.pi) * (
        np.cos(y_psr * 2.0 * np.pi) * np.sin(colatitude) * np.sin(theta_A) - np.cos(colatitude) * np.cos(theta_A)
    )

    phase_ensemble = phase - first_axis_phase + y0_det
    wrapped_phase = phase_ensemble - np.round(phase_ensemble)
    x = np.interp(wrapped_phase, x_det - np.round(np.min(x_det)) - 1.0, x_psr)
    y = np.interp(wrapped_phase, y_det - np.round(np.min(y_det)) - 1.0, y_psr)

    psi1 = psi_first(y)
    psi2 = psi_second(x)

    cos_psi_vel_1 = np.sin(theta_A) * -np.cos(theta_N + y * 2.0 * np.pi) * np.sin(theta_M) + np.cos(theta_M) * np.cos(theta_A)
    cos_psi_vel_2 = np.sin(theta_A) * -np.cos(theta_N + x * 2.0 * np.pi) * np.sin(theta_M) + np.cos(theta_M) * np.cos(theta_A)
    epsilon1 = 1.0 / pair.gamma_lor / (1.0 - beta_bulk * cos_psi_vel_1)
    epsilon2 = 1.0 / pair.gamma_lor / (1.0 + beta_bulk * cos_psi_vel_2)

    counts1 = shape.count_rate(psi1, epsilon1)
    counts2 = shape.count_rate(psi2, epsilon2)

    return BeamPairSimulation(
        colatitude_rad=colatitude,
        radius_over_rlc=radius,
        velocity_over_c=beta_bulk,
        psi_first_deg=psi1,
        psi_second_deg=psi2,
        epsilon_first=epsilon1,
        epsilon_second=epsilon2,
        counts_first=counts1,
        counts_second=counts2,
    )


def simulate_four_beam_band(
    phase: Array,
    geometry: FourBeamGeometry,
    band_model: EnergyBandModel,
    interpolation_grid_size: int = 2200,
) -> FourBeamSimulation:
    """Simulate a four-beam phaseogram for one energy band."""
    phase = _as_float_array(phase)
    lower = simulate_beam_pair(phase, geometry.theta_A_deg, geometry.lower, band_model.lower_shape, interpolation_grid_size)
    higher = simulate_beam_pair(phase, geometry.theta_A_deg, geometry.higher, band_model.higher_shape, interpolation_grid_size)
    background = float(band_model.background)
    model = background + lower.counts_first + lower.counts_second + higher.counts_first + higher.counts_second
    return FourBeamSimulation(
        phase=phase,
        background=background,
        beam1=lower.counts_first,
        beam2=lower.counts_second,
        beam3=higher.counts_first,
        beam4=higher.counts_second,
        model=model,
        lower_pair=lower,
        higher_pair=higher,
    )


def geometry_from_flat_parameters(params: Mapping[str, float]) -> FourBeamGeometry:
    """Build :class:`FourBeamGeometry` from flat parameters.

    Accepted keys follow the cleaned notation but common legacy aliases are also
    supported, e.g. ``angleA_deg`` for ``theta_A_deg``.
    """

    def get(clean: str, legacy: str | None = None) -> float:
        if clean in params:
            return float(params[clean])
        if legacy is not None and legacy in params:
            return float(params[legacy])
        raise KeyError(f"Missing parameter: {clean}")

    theta_A = get("theta_A_deg", "angleA_deg")
    lower = BeamPairGeometry(
        r_sin_colat_over_rlc=get("d_sin_theta_B", "d1sin"),
        phase_separation=get("P2_minus_P1", "P2minusP1"),
        t0_phase=get("t0", "t___0"),
        gamma_lor=get("gamma_lor_1", "lor"),
        bulk_zenith_deg=get("theta_M_deg", "thetaM_deg"),
        bulk_azimuth_deg=get("theta_N_deg", "thetaN_deg"),
        beam_zenith_deg=get("theta_C_deg", "thetaC_deg"),
        beam_azimuth_deg=get("theta_Q_deg", "thetaQ_deg"),
    )
    higher = BeamPairGeometry(
        r_sin_colat_over_rlc=get("D_sin_phi_B", "D3SIN"),
        phase_separation=get("P4_minus_P3", "P4minusP3"),
        t0_phase=get("T0", "T___O"),
        gamma_lor=get("gamma_lor_3", "LOR"),
        bulk_zenith_deg=get("phi_M_deg", "phiM_deg"),
        bulk_azimuth_deg=get("phi_N_deg", "phiN_deg"),
        beam_zenith_deg=get("phi_C_deg", "phiC_deg"),
        beam_azimuth_deg=get("phi_Q_deg", "phiQ_deg"),
    )
    return FourBeamGeometry(theta_A_deg=theta_A, lower=lower, higher=higher)


def band_model_from_flat_parameters(params: Mapping[str, float], prefix: str) -> EnergyBandModel:
    """Build one band model from flat parameters.

    Parameters
    ----------
    params:
        Flat parameter mapping.
    prefix:
        Energy-band prefix such as ``BM``, ``LO``, ``HI`` or ``TP``.
    """
    p = prefix.rstrip("_")
    lower = BeamShape(
        log10_A=float(params[f"{p}_LogNorm0"]),
        eta=float(params[f"{p}_Eta1"]),
        psi_c_deg=float(params[f"{p}_Psi_c0"]),
        beta=float(params[f"{p}_Beta"]),
    )
    higher = BeamShape(
        log10_A=float(params[f"{p}_LogAmpO"]),
        eta=float(params[f"{p}_eTA1"]),
        psi_c_deg=float(params[f"{p}_pSI_wO"]),
        beta=float(params[f"{p}_bETA"]),
    )
    return EnergyBandModel(lower_shape=lower, higher_shape=higher, background=float(params[f"{p}_Const"]))


def update_legacy_aliases(params: MutableMapping[str, float]) -> MutableMapping[str, float]:
    """Add paper-style aliases to a legacy parameter dictionary in-place."""
    aliases = {
        "angleA_deg": "theta_A_deg",
        "d1sin": "d_sin_theta_B",
        "P2minusP1": "P2_minus_P1",
        "t___0": "t0",
        "lor": "gamma_lor_1",
        "D3SIN": "D_sin_phi_B",
        "P4minusP3": "P4_minus_P3",
        "T___O": "T0",
        "LOR": "gamma_lor_3",
    }
    for legacy, clean in aliases.items():
        if legacy in params and clean not in params:
            params[clean] = params[legacy]
    return params
