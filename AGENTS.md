# AGENTS.md

Guidance for AI coding agents working in this repository.

## Project overview

CluE Oxide (Cluster Evolution) is a Rust program that simulates central
(electron) spin decoherence using Yang and Liu's cluster correlation expansion
(CCE). It reads a molecular structure (PDB, GRO, or other supported formats),
builds a nuclear spin bath around the detected electron spin, partitions the
bath into clusters, and computes pulse-sequence signals (e.g. Hahn echo).

- Crate: `clue_oxide` (binary and library), Rust edition 2021, GPL-3.0.
- User documentation: `manual/CluE_Oxide.tex` (built PDF at
  `manual/CluE_Oxide.pdf`). User-visible changes go in `change_log.txt`.
- Python bindings (pyCluE) live in `pyclue/`, a separate crate built with
  maturin and PyO3.

## Repository layout

```
src/
  main.rs            CLI entry point: parse args -> Config -> clue_oxide::run
  lib.rs             `run(config)`: time axis, output dir, calculate signals, write CSVs
  config.rs          `Config` struct, defaults, enums (ClusterMethod, PulseSequence, ...)
  config/
    config_toml.rs   serde TOML schema + ALLOWED_KEYS list (unknown keys are errors)
    particle_config.rs  per-group / per-isotope particle settings
    command_line_input.rs  CLI flags (-h, -l, -V, -W, -H, -O)
  structure/         structure input (pdb.rs, gro.rs), particle filters,
                     exchange groups, extended (periodic) structures
  cluster/           adjacency lists, cluster finding, partitioning
                     (partition.rs, methyl_clusters.rs), cluster TOML I/O
  quantum/           spin Hamiltonian construction and coupling tensors
  signal/            signal calculation, CCE, analytic 2-cluster signals
  clue_errors.rs     `CluEError` enum: the single error type for the crate
  isotopes.rs, elements.rs, physical_constants.rs  physical data tables
  space_3d.rs, math.rs, kmeans.rs, integration_grid.rs, symmetric_list_2d.rs
  info/              --help, --license, --version, --warranty text
tests/test_1omp.rs   integration test: compares a full run to a reference signal
assets/              structures, TOML configs, and reference CSVs used by tests
examples/            example inputs and plotting scripts (not maintained by tests)
pyclue/src/          current Python bindings (py_*.rs wrappers around the core crate)
pyclue/dev_src/, pyclue/old_src/  older/experimental copies; do not edit unless asked
manual/              LaTeX manual
```

## Build and test

Run all commands from the repository root: unit tests and the integration
test load files via relative paths such as `assets/TEMPO.pdb`.

```sh
cargo build --release     # binary at target/release/clue_oxide
cargo test                # unit tests (in-module #[cfg(test)]) + tests/test_1omp.rs
./clean_test_output.sh    # removes CluE-* output directories created by runs/tests
```

- The dev profile uses `opt-level = 3`, so debug builds are slow to compile but
  run tests at near-release speed.
- Linear algebra uses `ndarray-linalg` with `openblas-system` on Unix; the
  system needs OpenBLAS (and LAPACK) installed. `Cargo.toml` contains a
  commented alternative for a static OpenBLAS build and a Windows/Intel MKL
  section.
- When no `output_directory` is set, a run writes to `CluE-<hash>/`; these
  directories are gitignored.

### Python bindings

```sh
cd pyclue
maturin build             # wheel in pyclue/target/wheels
```

The PyO3 module is named `clue_oxide` (`module-name = "pyclue.clue_oxide"`).
When changing a public API in the core crate, check `pyclue/src/` for callers;
pyclue is a separate crate and is not compiled by `cargo test` at the root.

## Input format

The CLI takes a single TOML config file (`clue_oxide config.toml`). See
`assets/1omp_K26R1_0.003228966616048703.toml` for a representative example. The
old `.clue` token format was removed in 0.2.3-alpha.1; files such as
`examples/*/sim.clue` use the old format and no longer run.

When adding or renaming a config key:
1. Add the field to the serde structs in `src/config/config_toml.rs`.
2. Add the key to `ALLOWED_KEYS` in the same file and update the array length
   (`[&str; N]`). Unknown keys are reported as errors, so a missing entry
   breaks valid configs.
3. Handle it in `Config::set_from_config_toml` / `set_defaults` in `src/config.rs`.
4. Document it in `manual/CluE_Oxide.tex` and `change_log.txt`.

## Code conventions

Match the existing style rather than running `rustfmt` over whole files (the
code is not rustfmt-formatted and a reformat would produce large diffs):

- Two-space indentation; opening braces on the same line, often without a
  preceding space (`fn foo(){`, `match x{`).
- Explicit turbofish on collection types: `Vec::<f64>`, `Complex::<f64>`.
- Functions are separated by `//----------...` comment rules inside `impl`
  blocks; public items get a short `///` doc comment ("This function ...").
- Errors: return `Result<_, CluEError>` and propagate with `?`. To add an error,
  add a variant to `CluEError` in `src/clue_errors.rs` **and** a message arm in
  its `fmt` implementation.
- Randomness flows through a `ChaCha20Rng` seeded from `config.rng_seed`;
  pass the rng down rather than creating new generators so runs stay
  reproducible.
- Parallelism uses `rayon`.
- Unit tests live in `#[cfg(test)] mod tests` at the bottom of each module and
  frequently use files from `assets/`.
- Keep the code free of compiler warnings (recent commits cleaned these up).

## Versioning and releases

The version is defined in two places that must stay in sync:
- `Cargo.toml` (`version = ...`)
- `src/info/version.rs` (`print_version`)

Add an entry at the top of `change_log.txt` describing ADDED / CHANGED /
REMOVED items for user-visible changes.

## Numerical correctness

`tests/test_1omp.rs` compares a full simulation against
`assets/1omp_K26R1_0.003228966616048703_signal.csv` with an RMSD threshold of
`1e-12`. Changes to physics, constants, tensor construction, clustering, or
default settings can break this test; do not regenerate the reference data
unless the change in results is intended and explained.

