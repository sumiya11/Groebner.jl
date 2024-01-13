## Benchmarks

Benchmarks for `Groebner.jl`, `Singular`, `Maple`, `OpenF4`, and `msolve`.

The definitions of benchmark systems can be found in `benchmark_systems.jl`.

Computed Groebner bases will be verified against short certificates that are
assumed to be correct.

For running the benchmarks, you will need a Julia client v1.6+ installed.
See the official installation guide at
https://julialang.org/downloads/platform/.

## Running benchmarks

To benchmark all available software at once, execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl ALL
```

Also see below for software-specific instructions.

It is possible to specify command-line options. 
These are available for all software benchmarked here.
For example, the following command

```
julia one_script_to_run_them_all.jl ALL --timeout=600 --nworkers=20 --nruns=3 --validate=yes --suite=3
```

runs benchmarks under the following conditions:

- Timeout `600 s` for each benchmark.
- Distribute benchmarks over `20` worker processes.
- Re-run each benchmark `3` times and aggregate timings.
- Validate the correctness of resulting Groebner bases (may be slower).
- Use the `3`-rd benchmark suite (see `benchmark_systems.jl` for details).

## Groebner.jl

#### To run only `Groebner.jl` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl groebner
```

We compute the bases using the function `Groebner.groebner` with the option `threaded=:no`.

## Singular

*Running Singular benchmarks on Windows is not possible.*

#### To run only `Singular` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl singular
```

## Maple

#### For `Maple` benchmarks, you will need:

1. A Maple client installed. See the official installation guide at https://www.maplesoft.com/support/install/2021/Maple/Install.html

#### To run only `Maple` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl maple --bmaple=/path/to/maple/binary
```

Internally, this will use the `Groebner[Basis]` command with `method=fgb`.

## OpenF4

#### For `OpenF4` benchmarks, you will need:

1. OpenF4 lib installed. See the official installation guide at http://nauotit.github.io/openf4

#### To run only `OpenF4` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl openf4 --lopenf4=/path/to/openf4/lib
```

where `/path/to/openf4/lib` is the path to the directory where openf4 library is installed.

## msolve

#### For `msolve` benchmarks, you will need:

1. `msolve` installed. See the official installation guide at https://msolve.lip6.fr/ and https://github.com/algebraic-solving/msolve/blob/master/INSTALL

#### To run only `msolve` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl msolve --bmsolve=/path/to/msolve/binary
```

Internally, this will run `msolve -g 2 -l 44 -c 0 -f system.in -o /dev/null` for each benchmark system over finite fields and `msolve -g 2 -c 0 -f system.in -o /dev/null` for each benchmark system over the rationals.
