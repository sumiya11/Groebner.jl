
## Benchmarks for `Groebner.jl`

This directory contains benchmark systems and scripts to run them in `Groebner.jl`, `Singular`, and `Maple`.

All of the benchmarks are reproducible; the instructions to run benchmarks are listed in the corresponding directories.

- The `github.readme` directory is used to produce the table from the project Github readme https://github.com/sumiya11/Groebner.jl 

- The `paper` directory is used to produce tables from the paper "..."

- The `experiments` directory is used for development purposes. Benchmarks there are not reproducible, and the directory is a mess. 

- The `systems` directory contains test systems files. These include: 

    - `systems/biomodels`, polynomial chemical reaction network models obtained from https://odebase.org/

    - `systems/standard`, a short list of mostly zero-dimensional systems obtained from various sources

    - `systems/MQ`, a set of MQ problems obtained from https://www.mqchallenge.org/

    - `systems/for_gleb`, a couple of systems used in structural identifiability problems (see https://github.com/SciML/StructuralIdentifiability.jl/ for details)


