## Benchmarks

This directory can be used to run benchmarks for `Groebner.jl`, `Singular`, `Maple`, `OpenF4`, and `msolve`.

Benchmark systems are defined in `generate/benchmark_systems.jl`.

Benchmark results will be printed to the stdout and also written to a separate table in the `results` directory.

## Groebner.jl

#### For `Groebner.jl` benchmarks, you will need:

1. A Julia client v1.6+ installed. See the official installation guide at https://julialang.org/downloads/platform/


#### To run `Groebner.jl` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
julia one_script_to_rule_them_all.jl groebner
```

## Singular

#### For `Singular` benchmarks, you will need:

1. A Julia client v1.6+ installed. See the official installation guide at https://julialang.org/downloads/platform/

*Running Singular benchmarks on Windows platforms is not possible.*

#### To run `Singular` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
julia one_script_to_rule_them_all.jl singular
```

## Maple

#### For `Maple` benchmarks, you will need:

1. A Maple client v2021+ installed and added to `PATH`. See the official installation guide at https://www.maplesoft.com/support/install/2021/Maple/Install.html

#### To run `Maple` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
julia one_script_to_rule_them_all.jl maple
```

Internally, this will 
- Run `maple system.mpl` for each benchmark system
- TODO: the command

## OpenF4

#### For `OpenF4` benchmarks, you will need:

1. OpenF4 installed. See the official installation guide at http://nauotit.github.io/openf4

#### To run `OpenF4` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
julia one_script_to_rule_them_all.jl openf4 path_to_openf4_lib
```

Internally, this will 
- Run `g++ system.cpp -o bench & ./bench` for each benchmark system. TODO

## msolve

#### For `msolve` benchmarks, you will need:

1. `msolve` installed and added to `PATH`. See the official installation guide at https://github.com/algebraic-solving/msolve/blob/master/INSTALL

#### To run `msolve` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
julia one_script_to_rule_them_all.jl msolve
```

Internally, this will 
- Run `msolve -g system.in system.out` for each benchmark system.
