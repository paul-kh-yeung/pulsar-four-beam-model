# pulsar-four-beam-model

Reference implementation of the four-beam, geometry-first model used in

> **A common four-beam geometry reveals altitude-stratified GeV pulses in canonical young pulsars**  
> Paul K. H. Yeung and Takayuki Saito, accepted by *The Astrophysical Journal*.

The model decomposes Fermi-LAT gamma-ray phaseograms into two antipodal beam pairs:

- **Beams 1 and 2**: lower-altitude pair, close to the light-cylinder boundary;
- **Beams 3 and 4**: higher-altitude pair, contributing bridge/shoulder and ripple-like emission.

The code is intentionally mechanism-agnostic. It parametrizes emission-site locations,
bulk-flow directions, beam-axis orientations, phase-dependent Doppler shifts, and a
compact generalized-normal-distribution beam kernel.

## Status

This repository is intended as a clean public release of the analysis code. The importable 
implementation lives in `src/pulsar_four_beam_model/`.

The example Vela phaseogram files used by the default configuration `configs/vela_v1_97.json` are 
placed under `data/`.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## Quick start

After placing the Vela phaseogram files under `data/`, run:

```bash
python scripts/evaluate_solution.py configs/vela_v1_97.json
python scripts/plot_best_fit.py configs/vela_v1_97.json --output figures/vela_best_fit.png
python scripts/plot_psi_epsilon.py configs/vela_v1_97.json --output figures/vela_psi_epsilon.png

# Optional: rerun a compact Minuit fit.
# Fixed/free parameters are controlled by the JSON config; see the next section.
python scripts/fit_minuit.py configs/vela_v1_97.json --output results/vela_refit.json
```

The default Vela config expects these files:

```text
data/17yrVelaFermi_Data_0.06-0.15GeV_500.dat
data/17yrVelaFermi_Data_0.15-0.36GeV_500.dat
data/17yrVelaFermi_Data_0.36-0.9GeV_500.dat
data/17yrVelaFermi_Data_0.9-3GeV_500.dat
```

Each file should contain columns:

```text
phase  counts  [count_error]
```

If `count_error` is absent, Poisson errors `sqrt(counts)` are used.

## Data provenance

The Vela phaseogram files expected by the default configuration are derived from
Fermi-LAT observations. In the associated ApJ analysis, we used Pass 8 LAT
SOURCE-class events accumulated over 17 years of observations within a 3-degree
region of interest centred on the Vela pulsar, applying the standard good-time
selection and a zenith-angle cut to suppress Earth-limb contamination.

Photon arrival times were transformed to the Solar-system barycentre, and pulse
phases were assigned with the PINT timing package using an optimised Vela timing
ephemeris provided by Matthew Kerr. The resulting event list was then folded into
the phaseogram files used by this package.

The `.dat` files distributed or used with this repository are therefore not raw
Fermi-LAT event files. They are binned phaseograms derived from Fermi-LAT data
after barycentring and pulsar-phase assignment.

Please cite the associated ApJ paper when using these phaseograms or the model
implementation.

For full details on the provenance, timing, and phase-assignment procedure of the
example Vela phaseograms, see [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md).

## Controlling fixed and free parameters

`scripts/fit_minuit.py` reads the fixed/free parameter policy from the JSON config file.
The relevant keys in `configs/vela_v1_97.json` are:

```json
{
  "default_fixed": false,
  "fixed_parameters": [
    "BM_Alpha1",
    "BM_Alpha2",
    "BM_Eta2"
  ],
  "free_parameters": [],
  "strict_parameter_control": false
}
```

This mirrors the legacy Minuit syntax:

```python
m.fixed = False
m.fixed["BM_Alpha1"] = True
m.fixed["BM_Alpha2"] = True
m.fixed["BM_Eta2"] = True
```

The default policy follows the original research script: start with all active parameters
free (`default_fixed: false`), then fix the high-order/unused PLSEC coefficients listed
in `fixed_parameters`.

In this public package, the benchmark model is the compact 53-parameter generalized-normal-distribution
(GND) model used in the accepted paper. Therefore, many legacy PLSEC coefficients are
retained in the JSON file for provenance but are not active Minuit parameters in the
compact implementation. By default, `scripts/fit_minuit.py` skips such inactive names and
prints a short summary. Set

```json
"strict_parameter_control": true
```

if you want an inactive or misspelled parameter name to raise an error instead.

### Fix an active parameter

To fix an active parameter, add its name to `fixed_parameters`:

```json
"fixed_parameters": [
  "angleA_deg",
  "d1sin",
  "P2minusP1",
  "t___0"
]
```

### Force a parameter to be free

To force a parameter to be free, add it to `free_parameters`. This list is applied last,
so it can override a broader fixed list:

```json
{
  "default_fixed": true,
  "fixed_parameters": [],
  "free_parameters": [
    "BM_Const",
    "LO_Const",
    "HI_Const",
    "TP_Const"
  ]
}
```

This is useful for quick background-only refits.

### Active parameters in the compact 53-parameter model

The active parameters are the 17 global geometry/kinematic parameters plus 9 band-dependent
parameters for each energy band (`BM`, `LO`, `HI`, `TP`). These are the names that can be
meaningfully fixed or freed in the compact public model.

Global geometry/kinematic parameters:

```text
angleA_deg, d1sin, P2minusP1, t___0, lor,
thetaM_deg, thetaN_deg, thetaC_deg, thetaQ_deg,
D3SIN, P4minusP3, T___O, LOR,
phiM_deg, phiN_deg, phiC_deg, phiQ_deg
```

Band-dependent parameters, repeated for each band prefix `BM`, `LO`, `HI`, `TP`:

```text
{band}_LogNorm0, {band}_Eta1, {band}_Psi_c0, {band}_Beta,
{band}_LogAmpO, {band}_eTA1, {band}_pSI_wO, {band}_bETA,
{band}_Const
```

For example, the Vela config contains `BM_LogNorm0`, `BM_Eta1`, `BM_Psi_c0`,
`BM_Beta`, ..., `TP_Const`.

### Common fitting policies

Fix all global geometry/kinematic parameters, while refitting only the band-dependent
beam-shape and background parameters:

```json
"fixed_parameters": [
  "angleA_deg", "d1sin", "P2minusP1", "t___0", "lor",
  "thetaM_deg", "thetaN_deg", "thetaC_deg", "thetaQ_deg",
  "D3SIN", "P4minusP3", "T___O", "LOR",
  "phiM_deg", "phiN_deg", "phiC_deg", "phiQ_deg"
]
```

Fix everything except the four background constants:

```json
{
  "default_fixed": true,
  "fixed_parameters": [],
  "free_parameters": ["BM_Const", "LO_Const", "HI_Const", "TP_Const"]
}
```

Fix everything except the four beam-pair amplitudes and Doppler-response indices:

```json
{
  "default_fixed": true,
  "fixed_parameters": [],
  "free_parameters": [
    "BM_LogNorm0", "BM_Eta1", "BM_LogAmpO", "BM_eTA1",
    "LO_LogNorm0", "LO_Eta1", "LO_LogAmpO", "LO_eTA1",
    "HI_LogNorm0", "HI_Eta1", "HI_LogAmpO", "HI_eTA1",
    "TP_LogNorm0", "TP_Eta1", "TP_LogAmpO", "TP_eTA1"
  ]
}
```

After a fit, `scripts/fit_minuit.py` writes a `fit_control_summary` block to the output
JSON file, recording which active parameters were fixed/free and which legacy parameter
names were skipped.

## Core model notation

The notation follows the accepted ApJ paper as closely as possible. Some config names are
kept close to the legacy script for provenance.

| Config / code name | Paper notation | Meaning |
| --- | --- | --- |
| `angleA_deg` | `Θ_A` | angle between spin axis and line of sight |
| `d1sin` | `d sin θ_B / R_LC` | lower-altitude cylindrical radius |
| `D3SIN` | `D sin φ_B / R_LC` | higher-altitude cylindrical radius |
| `t___0`, `T___O` | `t_0`, `T_0` | detector-frame phases when sites pass closest to the line of sight |
| `lor`, `LOR` | `γ_Lor,1`, `γ_Lor,3` | effective bulk Lorentz factors |
| `thetaM_deg`, `thetaN_deg` | `θ_M`, `θ_N` | lower-pair bulk-flow direction |
| `phiM_deg`, `phiN_deg` | `φ_M`, `φ_N` | higher-pair bulk-flow direction |
| `thetaC_deg`, `thetaQ_deg` | `θ_C`, `θ_Q` | lower-pair beam-axis direction |
| `phiC_deg`, `phiQ_deg` | `φ_C`, `φ_Q` | higher-pair beam-axis direction |
| `{band}_LogNorm0`, `{band}_Eta1`, `{band}_Psi_c0`, `{band}_Beta` | `log10 A_1`, `η_1`, `Ψ_c,1`, `β_1` | band-dependent parameters for Beams 1 and 2 |
| `{band}_LogAmpO`, `{band}_eTA1`, `{band}_pSI_wO`, `{band}_bETA` | `log10 A_3`, `η_3`, `Ψ_c,3`, `β_3` | band-dependent parameters for Beams 3 and 4 |
| `{band}_Const` | `C_bkg` | constant background level in the phaseogram |

The detector-frame count-rate contribution from beam *i* is

```text
C_i(psi, epsilon) = A_i * epsilon**eta_i * exp[-(psi / Psi_c_i)**beta_i]
```

where `epsilon = E_bulk / E_det` is the inverse Doppler factor.

## Derived geometric quantities

The original legacy script also reports several additional geometric quantities derived
from the fitted parameters. These values are not independent fit parameters; they are
post-fit diagnostics useful for paper tables, sanity checks, and comparisons between the
lower-altitude and higher-altitude beam pairs.

| Quantity | Meaning |
| --- | --- |
| `d1_RLC` | spherical distance of the lower-altitude emission site from the pulsar centre, in units of `R_LC` |
| `thetaB_deg` | representative zenith angle of the lower-altitude emission site |
| `d1cos_RLC` | vertical coordinate of the lower-altitude emission site, `d1 cos(thetaB) / R_LC` |
| `P1phase`, `P2phase` | detector-frame closest-approach phases of Beams 1 and 2 |
| `D3_RLC` | spherical distance of the higher-altitude emission site from the pulsar centre, in units of `R_LC` |
| `phiB_deg` | representative zenith angle of the higher-altitude emission site |
| `D3cos_RLC` | vertical coordinate of the higher-altitude emission site, `D3 cos(phiB) / R_LC` |
| `P3phase`, `P4phase` | detector-frame closest-approach phases of Beams 3 and 4 |
| `Delta_theta_deg` | angular separation between the lower-altitude bulk-flow direction `(thetaM_deg, thetaN_deg)` and beam-axis direction `(thetaC_deg, thetaQ_deg)` |
| `Delta_phi_deg` | angular separation between the higher-altitude bulk-flow direction `(phiM_deg, phiN_deg)` and beam-axis direction `(phiC_deg, phiQ_deg)` |
| `AzDiff_site_deg` | azimuthal separation between the lower-altitude and higher-altitude emission sites |
| `AzDiff_beam_deg` | azimuthal separation between the corresponding lower-altitude and higher-altitude beam axes |

`Delta_theta_deg` and `Delta_phi_deg` quantify the small misalignment between
the fitted bulk-flow direction and the fitted beam-axis direction. They correspond
to the additional geometrical properties `Delta theta` and `Delta phi` reported in
the paper table. Both are computed using the scalar-product formula for two
directions written in spherical coordinates.

The legacy expressions are:

```python
P1phase = (
    t___0
    - thetaQ_deg / 360
    - d1sin / (2 * pi) * (cos(-thetaQ_deg * pi / 180) - 1) * sin(angleA_deg * pi / 180)
)
P2phase = P1phase + P2minusP1

P3phase = (
    T___O
    - phiQ_deg / 360
    - D3SIN / (2 * pi) * (cos(-phiQ_deg * pi / 180) - 1) * sin(angleA_deg * pi / 180)
)
P4phase = P3phase + P4minusP3

Delta_theta_deg = 180 / pi * acos(
    cos(thetaM_deg * pi / 180) * cos(thetaC_deg * pi / 180)
    + sin(thetaM_deg * pi / 180) * sin(thetaC_deg * pi / 180)
    * cos((thetaN_deg - thetaQ_deg) * pi / 180)
)
Delta_phi_deg = 180 / pi * acos(
    cos(phiM_deg * pi / 180) * cos(phiC_deg * pi / 180)
    + sin(phiM_deg * pi / 180) * sin(phiC_deg * pi / 180)
    * cos((phiN_deg - phiQ_deg) * pi / 180)
)

AzDiff_site_deg = 360 * (
    T___O
    - t___0
    + (P4minusP3 - P2minusP1) / 2
    + sin(angleA_deg * pi / 180) * (D3SIN - d1sin) / (2 * pi)
)
AzDiff_beam_deg = AzDiff_site_deg - (phiQ_deg - thetaQ_deg)
```

The recommended public-package implementation is to treat these as derived diagnostics,
not as part of the core likelihood. A clean place for the helper is:

```text
src/pulsar_four_beam_model/derived.py
```

and the natural scripts to report them are:

```text
scripts/evaluate_solution.py
scripts/fit_minuit.py
```

For reproducibility, a fit-output JSON can store them under a block such as:

```json
"derived_geometry": {
  "d1_RLC": 1.17,
  "thetaB_deg": 56.7,
  "d1cos_RLC": 0.64,
  "P1phase": 0.116,
  "P2phase": 0.548,
  "D3_RLC": 1.50,
  "phiB_deg": 66.3,
  "D3cos_RLC": 0.60,
  "P3phase": 0.220,
  "P4phase": 0.470,
  "Delta_theta_deg": 23.9,
  "Delta_phi_deg": 2.6,
  "AzDiff_site_deg": -12.4,
  "AzDiff_beam_deg": 38.5
}
```

## Repository layout

```text
src/pulsar_four_beam_model/
  model.py        # geometry, Doppler factor, four-beam simulation
  data.py         # phaseogram loading and weighted-Poisson helpers
  likelihood.py   # joint -2 log L
  plotting.py     # standard decomposition plots
  fit.py          # compact 53-parameter Minuit wrapper
  config.py       # JSON helpers
configs/
  vela_v1_97.json # Vela parameter set extracted from the legacy script
scripts/
  evaluate_solution.py
  plot_best_fit.py
  plot_psi_epsilon.py
  fit_minuit.py
docs/
  notation_mapping.md
DATA_PROVENANCE.md
ACKNOWLEDGMENTS.md
CITATION.cff
LICENSE
pyproject.toml
```

## Further development

Researchers interested in jointly developing, extending, or improving this
framework or its underlying algorithm are welcome to contact Paul K. H. Yeung
(<pkh91yg@icrr.u-tokyo.ac.jp>).

## Acknowledgements

Paul K. H. Yeung thanks Takayuki Saito for supervising this project. We are grateful to Dmitry Khangulyan for insightful discussions on phase-dependent
Doppler shifts. Sincere gratitude is given to Matthew Kerr for providing an
optimised timing ephemeris for the Vela pulsar. We thank Ievgen Vovk for useful
discussion.

## Citation

Please cite the associated ApJ paper and the archived Zenodo release of this
software when using this package. Citation metadata are provided in
[`CITATION.cff`](CITATION.cff).

## License

This project is distributed under the MIT License. See [`LICENSE`](LICENSE).
