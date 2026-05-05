"""Derived geometric quantities for the four-beam model.

These quantities are computed from the fitted parameters and are useful
for diagnostics, tables, and comparisons with the paper notation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import acos, atan, cos, floor, pi, sin
from typing import Mapping


@dataclass(frozen=True)
class DerivedGeometry:
    """Additional geometric quantities derived from fit parameters."""

    d1_RLC: float
    thetaB_deg: float
    d1cos_RLC: float
    P1phase: float
    P2phase: float

    D3_RLC: float
    phiB_deg: float
    D3cos_RLC: float
    P3phase: float
    P4phase: float

    Delta_theta_deg: float
    Delta_phi_deg: float
    AzDiff_site_deg: float
    AzDiff_beam_deg: float

    def to_dict(self) -> dict[str, float]:
        """Return the derived quantities as a plain dictionary."""
        return asdict(self)


def _angular_separation_deg(theta1_deg: float, phi1_deg: float,
                            theta2_deg: float, phi2_deg: float) -> float:
    """Angle between two directions given in spherical coordinates.

    Angles are in degrees.
    theta = zenith/polar angle from the spin axis
    phi   = azimuthal angle in the meridional-plane convention
    """
    theta1 = theta1_deg * pi / 180.0
    phi1 = phi1_deg * pi / 180.0
    theta2 = theta2_deg * pi / 180.0
    phi2 = phi2_deg * pi / 180.0

    cos_delta = (
        cos(theta1) * cos(theta2)
        + sin(theta1) * sin(theta2) * cos(phi1 - phi2)
    )

    # Numerical safety
    cos_delta = max(-1.0, min(1.0, cos_delta))

    return acos(cos_delta) * 180.0 / pi


def _pair_zenith_angle(radial_coordinate: float, phase_interval: float, angleA_rad: float) -> float:
    """Compute the representative zenith angle of one antipodal beam pair.

    This follows the legacy script convention used to recover thetaB / phiB
    from the fitted cylindrical radius and the observed phase interval.
    """
    wrapped_interval = phase_interval - floor(phase_interval)
    denominator = wrapped_interval - 0.5

    theta = atan(
        -radial_coordinate / (2.0 * pi) / denominator * 2.0 * cos(angleA_rad)
    )

    if theta < 0:
        theta += pi

    return theta


def compute_derived_geometry(parameters: Mapping[str, float]) -> DerivedGeometry:
    """Compute additional geometric quantities from a parameter dictionary.

    Expected parameter names follow the legacy/config convention:

    angleA_deg, d1sin, P2minusP1, t___0, thetaQ_deg,
    D3SIN, P4minusP3, T___O, phiQ_deg.
    """
    angleA_deg = float(parameters["angleA_deg"])
    d1sin = float(parameters["d1sin"])
    P2minusP1 = float(parameters["P2minusP1"])
    t___0 = float(parameters["t___0"])

    thetaM_deg = float(parameters["thetaM_deg"])
    thetaN_deg = float(parameters["thetaN_deg"])
    thetaC_deg = float(parameters["thetaC_deg"])
    thetaQ_deg = float(parameters["thetaQ_deg"])

    D3SIN = float(parameters["D3SIN"])
    P4minusP3 = float(parameters["P4minusP3"])
    T___O = float(parameters["T___O"])

    phiM_deg = float(parameters["phiM_deg"])
    phiN_deg = float(parameters["phiN_deg"])
    phiC_deg = float(parameters["phiC_deg"])
    phiQ_deg = float(parameters["phiQ_deg"])

    angleA_rad = angleA_deg * pi / 180.0
    thetaQ_rad = thetaQ_deg * pi / 180.0
    phiQ_rad = phiQ_deg * pi / 180.0

    P1phase = (
        t___0
        - thetaQ_deg / 360.0
        - d1sin / (2.0 * pi) * (cos(-thetaQ_rad) - 1.0) * sin(angleA_rad)
    )
    P2phase = P1phase + P2minusP1

    P3phase = (
        T___O
        - phiQ_deg / 360.0
        - D3SIN / (2.0 * pi) * (cos(-phiQ_rad) - 1.0) * sin(angleA_rad)
    )
    P4phase = P3phase + P4minusP3

    thetaB = _pair_zenith_angle(d1sin, P2phase - P1phase, angleA_rad)
    d1 = d1sin / sin(thetaB)

    phiB = _pair_zenith_angle(D3SIN, P4phase - P3phase, angleA_rad)
    D3 = D3SIN / sin(phiB)

    Delta_theta = _angular_separation_deg(
        thetaM_deg, thetaN_deg,
        thetaC_deg, thetaQ_deg,
    )

    Delta_phi = _angular_separation_deg(
        phiM_deg, phiN_deg,
        phiC_deg, phiQ_deg,
    )

    AzDiff_site = 360.0 * (
        T___O
        - t___0
        + (P4minusP3 - P2minusP1) / 2.0
        + sin(angleA_rad) * (D3SIN - d1sin) / (2.0 * pi)
    )

    AzDiff_beam = AzDiff_site - (phiQ_deg - thetaQ_deg)

    return DerivedGeometry(
        d1_RLC=d1,
        thetaB_deg=thetaB * 180.0 / pi,
        d1cos_RLC=d1 * cos(thetaB),
        P1phase=P1phase,
        P2phase=P2phase,
        D3_RLC=D3,
        phiB_deg=phiB * 180.0 / pi,
        D3cos_RLC=D3 * cos(phiB),
        P3phase=P3phase,
        P4phase=P4phase,
        Delta_theta_deg=Delta_theta,
        Delta_phi_deg=Delta_phi,
        AzDiff_site_deg=AzDiff_site,
        AzDiff_beam_deg=AzDiff_beam,
    )


def format_derived_geometry(derived: DerivedGeometry) -> str:
    """Format derived geometric quantities in the legacy print style."""
    return "\n".join(
        [
            f"d1 = {derived.d1_RLC:.8g} R_LC",
            f"thetaB = {derived.thetaB_deg:.8g} deg",
            f"d1cos = {derived.d1cos_RLC:.8g} R_LC",
            f"P1phase = {derived.P1phase:.8g}",
            f"P2phase = {derived.P2phase:.8g}",
            f"D3 = {derived.D3_RLC:.8g} R_LC",
            f"phiB = {derived.phiB_deg:.8g} deg",
            f"D3COS = {derived.D3cos_RLC:.8g} R_LC",
            f"P3phase = {derived.P3phase:.8g}",
            f"P4phase = {derived.P4phase:.8g}",
            f"Delta_theta = {derived.Delta_theta_deg:.8g} deg",
            f"Delta_phi = {derived.Delta_phi_deg:.8g} deg",
            f"AzDiff_site = {derived.AzDiff_site_deg:.8g} deg",
            f"AzDiff_beam = {derived.AzDiff_beam_deg:.8g} deg",
        ]
    )
