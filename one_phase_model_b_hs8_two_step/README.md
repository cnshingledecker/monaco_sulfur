# One-Phase Model B HS8 Two-Step Setup

This directory is a self-contained one-phase dense-cloud setup that keeps the validated `Model_B` physical conditions but replaces the sulfur chemistry with the reduced two-step HS8 hydrogenation network.

Purpose:

- preserve the exact one-phase `Model_B` conditions and initial abundances
- swap in the reduced two-step HS8 chemistry
- prepare a run-ready setup without executing the model yet

Model characteristics:

- `RHO = 1.0e4 cm^-3`
- `T = TDUST = 10 K`
- `AVST = AVIS = 10`
- `INIT_NON_ZERO = 1`
- `RADIOLYSIS = 1`
- `SUPRATHERMAL = 1`
- `FAST_BULK = 1`
- `SHINGLEDECKER_TUNN = 1`
- `TSTART = 1.0e2 yr`
- `TEND = 1.0e7 yr`

Chemistry changes relative to `one_phase_model_b_validation/`:

- `network.dat` includes `HS8` and `gHS8`
- `network.dat` adds `gH + gS8 -> gHS8`
- `network.dat` adds `gH + gHS8 -> gH2S + gS7`
- `network.dat` adds `HS8 <-> gHS8` freeze/desorb bookkeeping
- `model.inp` now tracks `HS8`, `gHS8`, and `bHS8` in the detailed-species list

Source-level requirement:

- the repository root `mod_calculate_rates.f90` now includes the fitted `SHINGLEDECKER_TUNN` special cases for `H + S8 -> HS8` and `H + HS8 -> H2S + S7` on both the surface and in the bulk

Files included:

- runtime inputs: `model.inp`, `network.dat`, `init_gas_ab.inp`, `init_surf_ab.inp`, `init_bulk_ab.inp`
- supporting data: `Lee_ea_17.txt`, `bulk_radiolysis.dat`, `class_2_suprathermal.dat`, `enthalpias.txt`, `radiolysis.dat`, `rd_eff.txt`

How to run later:

1. Build the executable in the repository root, for example with `make monaco`.
2. From this directory, run `./run_one_phase.sh`.

This directory intentionally omits generated outputs for now. It is a clean setup-only package.
