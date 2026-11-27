# Poly Mojo

Optional comparison backends for the [Poly Mojo](https://github.com/Spellbound-Mojo/poly_mojo)
polynomial library. They compute reduced monic Gröbner bases over QQ and native
prime fields, using degree reverse lexicographic order and the input variable order.
The bridge calls the public `polynomials.groebner` API directly.

## Prepare the environment

The [pinned Pixi environment](environment/pixi.toml) installs the published
community-channel **Poly Mojo 1.0.2** and **APN Mojo 1.0.2** packages with
**Mojo 1.2.0** and Python 3.12 on Linux x86-64. The lockfile records the exact
dependency builds. From the repository root:

```sh
pixi install --locked --manifest-path scripts/benchmark_different_software/generate/polymojo/environment/pixi.toml
eval "$(pixi shell-hook --manifest-path scripts/benchmark_different_software/generate/polymojo/environment/pixi.toml --shell bash)"
```

These dependencies are only needed when selecting a Poly Mojo backend. The
packages are installed precompiled; no library source build is needed.

During development, a source build can be supplied with `POLY_MOJO_PATH` pointing
to its precompiled package directory. The bridge never builds the library itself.
For an installed package, it finds the modules in the activated prefix's
`lib/mojo` directory automatically. `MOJO_COMMAND` can override the compiler on PATH;
retain the `MODULAR_HOME` supplied by environment activation.

From the repository root, build the bridge once:

```sh
export POLY_MOJO_DRIVER="$(python scripts/benchmark_different_software/generate/polymojo/build.py)"
```

The returned executable supports all variable counts, both coefficient fields,
and both certificate modes. The cache includes compiler, runtime, platform,
bridge source, and package hashes. No compilation occurs inside case timeouts.
Prime moduli and their nonnegative representatives must fit a signed 64-bit
integer; rational numerators and denominators have arbitrary precision.

The comparison runner's bundled Julia environments can be used as usual. To run
against the current Groebner.jl checkout instead, prepare a shared environment once:

```sh
JULIA_NUM_PRECOMPILE_TASKS=1 julia scripts/benchmark_different_software/generate/polymojo/setup.jl
export BENCHMARK_JULIA_ENV="$PWD/.local/environment"
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

This optional setup keeps package installation outside individual workers and
selects PrettyTables 2 for the runner's reporting API. Preserve the generated
Manifest with results. When `BENCHMARK_JULIA_ENV` is unset, existing environment
selection is unchanged.

## Run the upstream suites

```sh
julia scripts/benchmark_different_software/one_script_to_run_them_all.jl polymojo --suite=1 --nworkers=1 --nruns=3 --validate=yes
julia scripts/benchmark_different_software/one_script_to_run_them_all.jl polymojocertified --suite=1 --nworkers=1 --nruns=3 --validate=yes
```

Alternatively, pass `--bpolymojo=/absolute/path/to/driver` instead of setting
`POLY_MOJO_DRIVER`. Suite 1 uses GF(1073741827), suite 2 uses GF(1031), and suite 3
uses QQ. Suite 0 contains two small smoke cases. `ALL` attempts both optional
backends, skipping them through the existing setup-error handling when the
executable is unavailable.

## Timing and guarantees

| Backend | Timed computation | Separate verification |
|---|---|---|
| `polymojo` | `certify=False, trace=False` | None |
| `polymojocertified` | `certify=True, trace=True`, including retained proof evidence | Trace replay plus generator reduction and Buchberger's criterion |

Ordinary QQ reconstruction is probabilistic, like Groebner.jl's default
`certify=false`. Certified verification establishes ideal equality and the
Gröbner property. Full comparison with a reference reduced basis also checks
reducedness. The two production modes may choose different algorithms.

Each case performs one untimed warm-up and `--nruns` measured executions, reporting
the minimum computation time, as the Groebner.jl backend does. The certified
backend's standard report adds **Verify certificate, s**; the combined report adds
**polymojocertified (verify)**. Verification is excluded from computation time and
from fastest-computation highlighting. These are independent minima and should
not be added to claim a measured minimum total.

`samples.csv` retains each measured compute/verify pair, including full precision.
Input decoding and output serialization are untimed. The runner's process timeout
includes warm-up, every repetition, verification, and output. A result is reported
only after all requested repetitions and certificate checks succeed.

## Validate and test

Use the existing fingerprint checks for compatibility. For stronger validation,
compare complete outputs after running a reference backend with `--validate=yes`:

```sh
julia --project="$BENCHMARK_JULIA_ENV" scripts/benchmark_different_software/generate/polymojo/validate.jl REFERENCE_OUTPUT POLYMOJO_OUTPUT CERTIFIED_OUTPUT
```

This checks the coefficient field, ordered variable names, basis size, and every
coefficient and monomial, allowing only permutation of basis elements. It does
not rely on the runner's stored fingerprints. No fingerprint files are changed
by this comparison tool.

Some reference backends, including msolve over QQ, emit primitive integer
polynomials instead of monic polynomials. Use `validate.jl --normalize-monic
REFERENCE_OUTPUT POLYMOJO_OUTPUT` to allow nonzero scalar multiples of each basis
element. Field, variable order, basis size, and every normalized coefficient and
monomial must still agree. This option compares bases up to scalar normalization;
the default strict mode also checks that the output coefficients agree directly.

Tests against an activated Julia benchmark environment:

```sh
python -m unittest discover -s scripts/benchmark_different_software/generate/polymojo -p 'test_build.py'
julia --project="$BENCHMARK_JULIA_ENV" scripts/benchmark_different_software/generate/polymojo/test_validation.jl
julia --project="$BENCHMARK_JULIA_ENV" scripts/benchmark_different_software/generate/polymojo/test_adapter.jl
julia scripts/benchmark_different_software/one_script_to_run_them_all.jl polymojo --suite=0 --nworkers=1 --nruns=2 --validate=no
julia scripts/benchmark_different_software/one_script_to_run_them_all.jl polymojocertified --suite=0 --nworkers=1 --nruns=2 --validate=no
julia --project="$BENCHMARK_JULIA_ENV" scripts/benchmark_different_software/generate/polymojo/test_runner.jl
```

The setup tests use a fake compiler and need no Mojo installation. The remaining
tests exercise the real bridge, exact coefficients, different variable counts,
malformed input, optional-backend setup, and the standard reports.
