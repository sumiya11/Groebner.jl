
The directory is used to run benchmarks for `Groebner.jl`, `Singular`, and `Maple` for the arxiv preprint https://arxiv.org/abs/2304.06935.

*running OpenF4 benchmarks is not possible at the moment, work in progress*

#### For `Groebner.jl` and `Singular` benchmarks, you will need:

1. A Julia client v1.6+ installed. See the official installation guide at https://julialang.org/downloads/platform/


#### To run `Groebner.jl` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
> julia run_groebner.jl
```

#### To run `Singular` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
> julia run_singular.jl
```

*Running Singular benchmarks on Windows platforms is currently not possible.*

---

#### For `Maple` benchmarks, you will need:

1. A Maple client v2021+ installed. See the official installation guide at https://www.maplesoft.com/support/install/2021/Maple/Install.html


#### To run `Maple` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
> maple run_maple.mpl
```

