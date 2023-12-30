## Benchmarks

Benchmarks for `Groebner.jl`, `Singular`, `Maple`, `OpenF4`, and `msolve`. See
the instructions below for running benchmarks.

The definitions of benchmark systems can be found in `benchmark_systems.jl`.

Benchmark results will be printed to console and also written to a file in the
`results` directory.

Computed Groebner bases will be verified against short certificates that are
assumed to be correct.

## Groebner.jl

#### For `Groebner.jl` benchmarks, you will need:

1. A Julia client v1.6+ installed. See the official installation guide at https://julialang.org/downloads/platform/

#### To run `Groebner.jl` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl groebner
```

We run Groebner.jl benchmarks using the function `Groebner.groebner` with default options.

It is possible to specify command-line options. 
These are available for all benchmarked software. 
For example, the following command

```
julia one_script_to_run_them_all.jl groebner --timeout=600 --nworkers=20 --nruns=3 --validate=yes --benchmark=2
```

runs Groebner.jl benchmarks under the following conditions:

- Timeout `600 s` for each benchmark
- Distribute benchmarks over `20` worker processes
- Re-run each benchmark `3` times and aggregate timings
- Validate the correctness of resulting Groebner bases
- Use `2`-nd benchmark suite (see `benchmark_systems.jl` for details) 

## Singular

#### For `Singular` benchmarks, you will need:

1. A Julia client v1.6+ installed. See the official installation guide at https://julialang.org/downloads/platform/

*Running Singular benchmarks on Windows is not possible.*

#### To run `Singular` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl singular
```

## Maple

#### For `Maple` benchmarks, you will need:

1. A Maple client v2021+ installed and added to `PATH`. See the official installation guide at https://www.maplesoft.com/support/install/2021/Maple/Install.html

#### To run `Maple` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl maple
```

You can use the option `--bmaple` to specify the path to the Maple executable.

Internally, this will generate `system.mpl` and run `maple system.mpl` for each benchmark system.

## OpenF4

#### For `OpenF4` benchmarks, you will need:

1. OpenF4 installed. See the official installation guide at http://nauotit.github.io/openf4

#### To run `OpenF4` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl openf4 --lopenf4=/path/to/openf4/lib
```

where `/path/to/openf4/lib` is the path to the directory where openf4 library is installed.

Internally, this will generate `system.cpp`, compile `system.cpp` (with `g++`), and execute the resulting binary for each benchmark system.

## msolve

#### For `msolve` benchmarks, you will need:

1. `msolve` installed and added to `PATH`. See the official installation guide at https://msolve.lip6.fr/ and https://github.com/algebraic-solving/msolve/blob/master/INSTALL

#### To run `msolve` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl msolve
```

You can use the option `--bmsolve` to specify the path to the Maple executable.

Internally, this will run `msolve -g 2 -l 44 -c 0 -f system.in -o /dev/null` for each benchmark system.
