"""Plotting helpers for four-beam phaseogram decompositions."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .data import Phaseogram
from .model import EnergyBandModel, FourBeamGeometry, simulate_four_beam_band


def plot_phaseogram_decomposition(
    phaseogram: Phaseogram,
    geometry: FourBeamGeometry,
    band_model: EnergyBandModel,
    title: str = "",
    off_phase: Tuple[float, float] | None = (-0.2, 0.08),
    output: str | Path | None = None,
):
    """Plot observed counts and the four model components for one energy band."""
    sim = simulate_four_beam_band(phaseogram.phase, geometry, band_model)
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    if off_phase is not None:
        ax.fill_between(
            phaseogram.phase,
            0,
            1,
            where=(phaseogram.phase >= off_phase[0]) & (phaseogram.phase <= off_phase[1]),
            transform=ax.get_xaxis_transform(),
            color="0.9",
            zorder=0,
        )
    phases = geometry.axis_phases()
    for key, color in zip(["beam1", "beam2", "beam3", "beam4"], ["C0", "C1", "C4", "C3"]):
        ax.axvline(phases[key], ls=":", color=color, lw=1.0)
    ax.axhline(sim.background, ls=":", color="0.2", lw=1.0)
    ax.errorbar(phaseogram.phase, phaseogram.counts, phaseogram.count_errors, fmt=".", color="k", ms=3, lw=0.8, label="Data")
    ax.plot(phaseogram.phase, sim.model, color="0.25", lw=2.0, label="Total")
    ax.plot(phaseogram.phase, sim.background + sim.beam1, color="C0", lw=1.1, label="Bkg+Beam1")
    ax.plot(phaseogram.phase, sim.background + sim.beam2, color="C1", lw=1.1, label="Bkg+Beam2")
    ax.plot(phaseogram.phase, sim.background + sim.beam3, color="C4", lw=1.1, label="Bkg+Beam3")
    ax.plot(phaseogram.phase, sim.background + sim.beam4, color="C3", lw=1.1, label="Bkg+Beam4")
    ax.set_xlabel("Observed phase")
    ax.set_ylabel("Counts")
    if title:
        ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, dpi=200)
    return fig, ax


def plot_all_bands(
    phaseograms: Mapping[str, Phaseogram],
    geometry: FourBeamGeometry,
    band_models: Mapping[str, EnergyBandModel],
    band_order: Sequence[str],
    titles: Mapping[str, str] | None = None,
    off_phase: Tuple[float, float] | None = (-0.2, 0.08),
    output: str | Path | None = None,
):
    """Plot stacked four-beam decompositions for multiple energy bands."""
    n = len(band_order)
    fig, axes = plt.subplots(
        n,
        1,
        figsize=(5.0, 2.5 * n),
        sharex=True,
        gridspec_kw={
            "wspace": 0,
            "hspace": 0,
            "height_ratios": [1] * n,
        },
    )

    if n == 1:
        axes = [axes]

    titles = titles or {}
    phases = geometry.axis_phases()
    colors = {"beam1": "C0", "beam2": "C1", "beam3": "C4", "beam4": "C3"}

    for idx, (ax, band) in enumerate(zip(axes, band_order)):
        pg = phaseograms[band]
        sim = simulate_four_beam_band(pg.phase, geometry, band_models[band])

        if off_phase is not None:
            ax.fill_between(
                pg.phase,
                0,
                1,
                where=(pg.phase >= off_phase[0]) & (pg.phase <= off_phase[1]),
                transform=ax.get_xaxis_transform(),
                facecolor="gainsboro",
                zorder=0,
                label="OFF phase" if idx == 0 else None,
            )

        for key, color in colors.items():
            ax.axvline(phases[key], ls=":", color=color, lw=1.0)

        ax.axhline(
            sim.background,
            ls=":",
            color="0.2",
            lw=1.0,
            label="Bkg" if idx == 0 else None,
        )
        ax.plot(
            pg.phase,
            sim.model,
            color="0.25",
            lw=2.0,
            zorder=1,
            label="Total" if idx == 0 else None,
        )
        ax.errorbar(
            pg.phase,
            pg.counts,
            pg.count_errors,
            fmt=".",
            color="k",
            ms=3,
            lw=0.8,
            zorder=2,
            label="Data" if idx == 0 else None,
        )
        ax.plot(
            pg.phase,
            sim.background + sim.beam1,
            color="C0",
            lw=1.1,
            zorder=3,
            label="Bkg+Beam1" if idx == 0 else None,
        )
        ax.plot(
            pg.phase,
            sim.background + sim.beam2,
            color="C1",
            lw=1.1,
            zorder=4,
            label="Bkg+Beam2" if idx == 0 else None,
        )
        ax.plot(
            pg.phase,
            sim.background + sim.beam3,
            color="C4",
            lw=1.1,
            zorder=5,
            label="Bkg+Beam3" if idx == 0 else None,
        )
        ax.plot(
            pg.phase,
            sim.background + sim.beam4,
            color="C3",
            lw=1.1,
            zorder=6,
            label="Bkg+Beam4" if idx == 0 else None,
        )

        ax.set_ylabel(f"{titles.get(band, band)} Counts")

    axes[0].legend(fontsize=8, loc="best")
    axes[-1].set_xlabel("Observed phase")

    fig.tight_layout()

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=200)

    return fig, axes


def plot_psi_epsilon(
    phaseogram: Phaseogram,
    geometry: FourBeamGeometry,
    band_model: EnergyBandModel,
    title: str = "",
    off_phase: Tuple[float, float] | None = (-0.2, 0.08),
    output: str | Path | None = None,
):
    """Plot beam-axis angle psi and inverse Doppler factor epsilon.

    This reproduces the diagnostic plot used in the legacy script:

    - top panel: beam-axis--line-of-sight angle psi;
    - bottom panel: inverse Doppler factor epsilon = E_bulk / E_det.

    The quantities are determined by the fitted geometry and bulk-flow
    parameters. They are effectively independent of the beam-shape
    normalisations, but a band model is passed in so that the same simulation
    interface can be reused.
    """
    sim = simulate_four_beam_band(phaseogram.phase, geometry, band_model)
    phases = geometry.axis_phases()

    fig, ax = plt.subplots(
        2,
        sharex=True,
        figsize=(4.0, 5.0),
        gridspec_kw={"wspace": 0, "hspace": 0, "height_ratios": [1, 1]},
    )

    colors = {
        "beam1": "C0",
        "beam2": "C1",
        "beam3": "C4",
        "beam4": "C3",
    }

    for axis in ax:
        if off_phase is not None:
            axis.fill_between(
                phaseogram.phase,
                0,
                1,
                where=(phaseogram.phase >= off_phase[0])
                & (phaseogram.phase <= off_phase[1]),
                transform=axis.get_xaxis_transform(),
                facecolor="gainsboro",
                zorder=0,
            )

        for key, color in colors.items():
            axis.axvline(phases[key], ls=":", color=color, lw=1.0)

    ax[0].plot(
        phaseogram.phase,
        sim.lower_pair.psi_first_deg,
        color=colors["beam1"],
        label="Beam1",
    )
    ax[0].plot(
        phaseogram.phase,
        sim.lower_pair.psi_second_deg,
        color=colors["beam2"],
        label="Beam2",
    )
    ax[0].plot(
        phaseogram.phase,
        sim.higher_pair.psi_first_deg,
        color=colors["beam3"],
        label="Beam3",
    )
    ax[0].plot(
        phaseogram.phase,
        sim.higher_pair.psi_second_deg,
        color=colors["beam4"],
        label="Beam4",
    )
    ax[0].set_ylabel(r"Beam-axis--l.o.s. angle $\psi$ ($^\circ$)")
    ax[0].legend(fontsize=8, loc="best")

    ax[1].plot(
        phaseogram.phase,
        sim.lower_pair.epsilon_first,
        color=colors["beam1"],
        label="Beam1",
    )
    ax[1].plot(
        phaseogram.phase,
        sim.lower_pair.epsilon_second,
        color=colors["beam2"],
        label="Beam2",
    )
    ax[1].plot(
        phaseogram.phase,
        sim.higher_pair.epsilon_first,
        color=colors["beam3"],
        label="Beam3",
    )
    ax[1].plot(
        phaseogram.phase,
        sim.higher_pair.epsilon_second,
        color=colors["beam4"],
        label="Beam4",
    )
    ax[1].set_ylabel(r"$\epsilon = E_{\rm bulk}/E_{\rm det}$")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("Observed phase")

    if title:
        ax[0].set_title(title)

    fig.tight_layout()

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=200)

    return fig, ax
