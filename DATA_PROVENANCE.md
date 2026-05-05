# Data provenance

This document describes the provenance of the example Vela phaseogram files used
by the default configuration of `pulsar-four-beam-model`.

## Data product

The Vela `.dat` files expected by the default configuration are binned gamma-ray phaseograms
derived from Fermi Large Area Telescope (Fermi-LAT) observations of the Vela
pulsar, PSR J0835-4510.

They are **not** raw Fermi-LAT photon event files. They are post-processed,
phase-folded count histograms prepared for the four-beam modelling workflow.

The default configuration expects four Vela phaseogram files:

```text
data/17yrVelaFermi_Data_0.06-0.15GeV_500.dat
data/17yrVelaFermi_Data_0.15-0.36GeV_500.dat
data/17yrVelaFermi_Data_0.36-0.9GeV_500.dat
data/17yrVelaFermi_Data_0.9-3GeV_500.dat
```

These correspond to four detected-energy bands spanning 60 MeV to 3 GeV:

```text
0.06--0.15 GeV
0.15--0.36 GeV
0.36--0.9 GeV
0.9--3 GeV
```

Each file contains binned phaseogram data with columns:

```text
phase  counts  [count_error]
```

If the optional `count_error` column is absent, the package uses Poisson errors,
`sqrt(counts)`, when plotting or evaluating the phaseograms.

## Fermi-LAT event selection

For the associated ApJ analysis, the Vela phaseograms were derived from Pass 8
Fermi-LAT SOURCE-class events accumulated over 17 years of observations.

The event selection used a circular region of interest with radius 3 degrees
centred on the Vela pulsar. Good-time intervals were selected following the
standard Fermi-LAT recommendations, including a zenith-angle cut of less than
90 degrees to suppress contamination from the Earth's limb.

No per-photon probability weights were assigned in the four-beam modelling.
Unpulsed and background contributions are instead represented by a constant
background term fitted simultaneously with the pulsed beam components.

## Timing and phase assignment

Photon arrival times were transformed to the Solar-system barycentre.

Pulse phases were assigned with the PINT timing package using an optimised Vela
timing ephemeris provided by Matthew Kerr. The resulting phase-assigned event
list was then folded into the binned `.dat` phaseograms used by this package.

## Relation to this software package

When provided, the `.dat` files are inputs to the modelling scripts. 
They are not part of the core model implementation.

The model itself is implemented in:

```text
src/pulsar_four_beam_model/
```

The default Vela configuration is stored in:

```text
configs/vela_v1_97.json
```

Users who wish to reproduce the phaseograms from raw photon data should download
the relevant Fermi-LAT data from the Fermi Science Support Center and repeat the
barycentring and pulsar-phase assignment procedure using an appropriate timing
solution.

## Credit and citation

If you use these phaseograms or the associated modelling code, please cite the
associated ApJ paper:

> Paul K. H. Yeung and Takayuki Saito,
> "A common four-beam geometry reveals altitude-stratified GeV pulses in
> canonical young pulsars",
> The Astrophysical Journal, 2026.

Please also acknowledge the Fermi-LAT collaboration for operating the mission
and providing public access to the data.

For the Vela phaseograms specifically, please note that the optimised timing
ephemeris used for phase assignment was provided by Matthew Kerr.
