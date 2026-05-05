# Notation mapping from the legacy script to the paper-style code

The original script `UltimateToyModel_Quadruple_advanced_v1.97.py` was written as a
research notebook-like executable. The cleaned package keeps the same mathematical
structure but renames variables to match the ApJ paper.

## Global geometry

| Legacy name | Clean name | Paper notation | Description |
| --- | --- | --- | --- |
| `angleA_deg` | `theta_A_deg` | `Θ_A` | spin-axis inclination from the line of sight |
| `d1sin` | `d_sin_theta_B` | `d sin θ_B / R_LC` | lower-altitude cylindrical radius |
| `D3SIN` | `D_sin_phi_B` | `D sin φ_B / R_LC` | higher-altitude cylindrical radius |
| `P2minusP1` | `P2_minus_P1` | `P2 - P1` | lower-pair axis-phase separation |
| `P4minusP3` | `P4_minus_P3` | `P4 - P3` | higher-pair axis-phase separation |
| `t___0` | `t0` | `t_0` | phase at which lower-altitude site passes closest to line of sight |
| `T___O` | `T0` | `T_0` | phase at which higher-altitude site passes closest to line of sight |

## Lower-altitude beam pair: Beams 1 and 2

| Legacy name | Clean name | Paper notation |
| --- | --- | --- |
| `lor` | `gamma_lor_1` | `γ_Lor,1` |
| `thetaM_deg` | `theta_M_deg` | `θ_M` |
| `thetaN_deg` | `theta_N_deg` | `θ_N` |
| `thetaC_deg` | `theta_C_deg` | `θ_C` |
| `thetaQ_deg` | `theta_Q_deg` | `θ_Q` |

## Higher-altitude beam pair: Beams 3 and 4

| Legacy name | Clean name | Paper notation |
| --- | --- | --- |
| `LOR` | `gamma_lor_3` | `γ_Lor,3` |
| `phiM_deg` | `phi_M_deg` | `φ_M` |
| `phiN_deg` | `phi_N_deg` | `φ_N` |
| `phiC_deg` | `phi_C_deg` | `φ_C` |
| `phiQ_deg` | `phi_Q_deg` | `φ_Q` |

## Energy-band prefixes

| Prefix | Energy band in the default Vela config |
| --- | --- |
| `BM` | 0.06--0.15 GeV |
| `LO` | 0.15--0.36 GeV |
| `HI` | 0.36--0.9 GeV |
| `TP` | 0.9--3 GeV |

## Band-dependent beam-shape parameters

For each prefix `XX`:

| Legacy lower-pair parameter | Paper notation |
| --- | --- |
| `XX_LogNorm0` | `log10 A_1` |
| `XX_Eta1` | `η_1` |
| `XX_Psi_c0` | `Ψ_c,1` |
| `XX_Beta` | `β_1` |

| Legacy higher-pair parameter | Paper notation |
| --- | --- |
| `XX_LogAmpO` | `log10 A_3` |
| `XX_eTA1` | `η_3` |
| `XX_pSI_wO` | `Ψ_c,3` |
| `XX_bETA` | `β_3` |

The legacy parameters `Gamma*`, `Alpha*`, `Eta2`--`Eta4`, `gAMMA*`, `aLPHA*`, and
`eTA2`--`eTA4` are preserved in the config for provenance but are not used by the
clean benchmark implementation, because the accepted-paper model uses the compact
GND expression with those higher-order terms fixed to zero.
