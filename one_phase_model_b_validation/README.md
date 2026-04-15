# One-Phase Model B Validation

This directory is a self-contained one-phase dense-cloud sulfur validation case.

Purpose:

- reproduce the archived paper-era `Model_B` setup as a one-phase run
- confirm that the benchmark `bS8` abundance curve comes from the root `model.inp` plus root `init_*` files, not from the archived three-phase handoff chain

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

Files included:

- runtime inputs: `model.inp`, `network.dat`, `init_gas_ab.inp`, `init_surf_ab.inp`, `init_bulk_ab.inp`
- supporting data: `Lee_ea_17.txt`, `bulk_radiolysis.dat`, `class_2_suprathermal.dat`, `enthalpias.txt`, `radiolysis.dat`, `rd_eff.txt`
- validation target: `expected_outputs/bS8.ab`

How to run:

1. Build the executable in the repository root, for example with `make monaco`.
2. From this directory, run `./run_one_phase.sh`.

Expected validation result:

- the generated `ab/bS8.ab` should match `expected_outputs/bS8.ab` byte-for-byte
- the reproduced curve peaks at `1.41e-7` near `6.31e6 yr` and declines to `1.19e-11` at `1.0e7 yr`

This directory intentionally does not include the full generated output tree, to keep the branch lightweight and reviewable.
