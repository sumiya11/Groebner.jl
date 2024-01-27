## Benchmarks

Benchmarks for `Groebner.jl`, `Maple`, `msolve`, `OpenF4`, `Singular`.

The definitions of benchmark systems can be found in `benchmark_systems.jl`.

Computed bases will be verified against short SHA certificate, which are assumed to
be correct.

For running the benchmarks, you will need a Julia client v1.6+ installed (preferably 1.9).
See the official installation guide at
https://julialang.org/downloads/platform/.

## Running benchmarks

To benchmark all available software at once, execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl ALL
```

See below for software-specific instructions.

It is possible to specify command-line options. 
These options are available for all software benchmarked here.
For example, the following command

```
julia one_script_to_run_them_all.jl ALL --timeout=600 --nworkers=16 --nruns=3 --validate=no --suite=1 --bmaple=/path/to/maple/binary --lopenf4=/path/to/openf4/lib --bmsolve=/path/to/msolve/binary
```

runs benchmarks under the following conditions:

- Timeout `600 s` for each benchmark system.
- Distribute benchmarks over `16` worker processes.
- Re-run each benchmark `3` times and average the timings.
- Do not validate the correctness of resulting Groebner bases.
- Use the `1`-st benchmark suite (see `benchmark_systems.jl` for details).
- Use Maple binary from `/path/to/maple/binary`.
- Use OpenF4 library installation from `/path/to/openf4/lib/`.
- Use msolve binary from `/path/to/msolve/binary`.

You can still run the command even if you do not have some of these software. 
For example, if you do not have `msolve` and `OpenF4`, omit the `--bmsolve` and `--lopenf4` options, and the script will benchmark only `Groebner.jl`, `Maple`, and `Singular`.
We use `Singular.jl`, thus the only necessary thing to run Singular benchmarks is a Julia client.

It is advisable to run the benchmarks with no other maple/msolve tasks running in background in parallel. This is because the benchmark script eagerly cleans up the processes and may accidentally kill an innocent maple/msolve process.

## Groebner.jl

#### To run only `Groebner.jl` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl groebner
```

We compute the bases using the function `Groebner.groebner` with the option `threaded=:no` and other options with default values.

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

1. OpenF4 library installed. See the official installation guide at http://nauotit.github.io/openf4

#### To run only `OpenF4` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl openf4 --lopenf4=/path/to/openf4/lib
```

where `/path/to/openf4/lib` is the path to the directory where the `OpenF4` library is installed.

## msolve

#### For `msolve` benchmarks, you will need:

1. `msolve` client installed. See the official installation guide at https://msolve.lip6.fr/ and https://github.com/algebraic-solving/msolve/blob/master/INSTALL

#### To run only `msolve` benchmarks

1. Execute this command in your favorite terminal from this directory:

```
julia one_script_to_run_them_all.jl msolve --bmsolve=/path/to/msolve/binary
```

Internally, this will run

```
msolve -g 2 -l 44 -c 0 -f system.in -o /dev/null
```

for each benchmark system over finite fields, where `msolve` is equal to `/path/to/msolve/binary`.
For benchmarks over the rationals, we do not specify the linear algebra option `-l`.
