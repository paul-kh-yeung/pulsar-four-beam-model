#!/usr/bin/env python
"""Run a compact 53-parameter Minuit fit from a saved config.

This is intentionally a lightweight public wrapper. For publication-scale scans,
use this as a reproducible starting point and keep long random-restart campaigns
in a separate workflow script.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from iminuit import Minuit

from pulsar_four_beam_model.config import read_json, write_json
from pulsar_four_beam_model.data import load_phaseogram
from pulsar_four_beam_model.fit import build_minuit
from pulsar_four_beam_model.derived import (
    compute_derived_geometry,
    format_derived_geometry,
)


def _normalise_name_list(value: Any, key: str) -> list[str]:
    """Return a JSON value as a list of parameter names."""
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise TypeError(f"`{key}` must be a list of parameter-name strings.")
    return value


def apply_fixed_free_settings(minuit: Minuit, cfg: dict[str, Any]) -> dict[str, Any]:
    """Apply fixed/free controls from a JSON config to a Minuit object.

    The config keys mirror the legacy script idiom:

        m.fixed = False
        m.fixed["BM_Alpha1"] = True

    Equivalent JSON keys:

        "default_fixed": false,
        "fixed_parameters": ["BM_Alpha1", ...],
        "free_parameters": []

    Notes
    -----
    The public benchmark implementation is the compact 53-parameter GND model.
    Some legacy coefficients are retained in the JSON file for provenance, but
    are not active Minuit parameters in this compact model. By default, such
    inactive names are skipped and reported. Set "strict_parameter_control": true
    in the JSON file to make unknown names raise an error instead.
    """
    active_parameters = set(minuit.parameters)

    default_fixed = bool(cfg.get("default_fixed", False))
    fixed_parameters = _normalise_name_list(
        cfg.get("fixed_parameters", []),
        "fixed_parameters",
    )
    free_parameters = _normalise_name_list(
        cfg.get("free_parameters", []),
        "free_parameters",
    )
    strict = bool(cfg.get("strict_parameter_control", False))

    # Mirrors legacy `m.fixed = False` or `m.fixed = True`.
    minuit.fixed = default_fixed

    skipped_fixed: list[str] = []
    skipped_free: list[str] = []

    for name in fixed_parameters:
        if name in active_parameters:
            minuit.fixed[name] = True
        elif strict:
            raise KeyError(f"Parameter listed in fixed_parameters is not active: {name}")
        else:
            skipped_fixed.append(name)

    # Apply explicit frees last, so they can override a broad fixed list.
    for name in free_parameters:
        if name in active_parameters:
            minuit.fixed[name] = False
        elif strict:
            raise KeyError(f"Parameter listed in free_parameters is not active: {name}")
        else:
            skipped_free.append(name)

    active_fixed = [name for name in minuit.parameters if minuit.fixed[name]]
    active_free = [name for name in minuit.parameters if not minuit.fixed[name]]

    return {
        "default_fixed": default_fixed,
        "n_active_parameters": len(minuit.parameters),
        "n_active_fixed": len(active_fixed),
        "n_active_free": len(active_free),
        "active_fixed_parameters": active_fixed,
        "active_free_parameters": active_free,
        "skipped_fixed_parameters": skipped_fixed,
        "skipped_free_parameters": skipped_free,
        "strict_parameter_control": strict,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument("--data-root", default=".", help="Root directory prepended to data paths in the config")
    parser.add_argument("--output", default="fit_result.json", help="Output JSON file")
    parser.add_argument("--max-calls", type=int, default=0, help="Maximum Minuit calls; 0 lets Minuit choose")
    args = parser.parse_args()

    cfg = read_json(args.config)
    bands = cfg.get("band_order", list(cfg["bands"].keys()))
    data_root = Path(args.data_root)
    phaseograms = {
        band: load_phaseogram(data_root / cfg["bands"][band]["file"])
        for band in bands
    }
    minuit = build_minuit(phaseograms, cfg["parameters"], bands=bands)
    fit_control = apply_fixed_free_settings(minuit, cfg)

    print(
        "Parameter control: "
        f"{fit_control['n_active_free']} free / "
        f"{fit_control['n_active_fixed']} fixed "
        f"out of {fit_control['n_active_parameters']} active parameters."
    )

    if fit_control["skipped_fixed_parameters"] or fit_control["skipped_free_parameters"]:
        print(
            "Skipped inactive legacy parameter names: "
            f"{len(fit_control['skipped_fixed_parameters'])} fixed-list, "
            f"{len(fit_control['skipped_free_parameters'])} free-list."
        )

    if args.max_calls > 0:
        minuit.migrad(ncall=args.max_calls)
    else:
        minuit.migrad()
    try:
        minuit.hesse()
    except Exception as exc:  # noqa: BLE001
        print(f"HESSE failed: {exc}")

    result = dict(cfg)
    result["fit_status"] = {
        "fval": float(minuit.fval),
        "valid": bool(minuit.valid),
        "accurate": bool(minuit.accurate),
    }
    result["fit_control_summary"] = fit_control
    result["parameters"] = dict(cfg["parameters"])
    result["parameters"].update(minuit.values.to_dict())
    derived = compute_derived_geometry(result["parameters"])

    print("\nDerived geometric quantities")
    print("----------------------------")
    print(format_derived_geometry(derived))

    result["derived_geometry"] = derived.to_dict()
    result["errors"] = minuit.errors.to_dict()
    write_json(result, args.output)
    print(f"Saved {args.output}")
    print(f"fval = {minuit.fval:.8g}; valid = {minuit.valid}")


if __name__ == "__main__":
    main()
