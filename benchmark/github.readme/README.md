
The directory used to run a tiny set of benchmark systems with `Groebner.jl`, `Singular`, and `Maple`.

#### For `Groebner.jl` and `Singular` benchmarks, you will need:

1. A Julia client v1.6+ installed. See the official installation guide at

https://julialang.org/downloads/platform/


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

*Running singular.jl benchmark on Windows platforms is currently not possible.*

---

#### For `Maple` benchmarks, you will need:

1. A Maple client v2021+ installed. See the official installation guide at

https://www.maplesoft.com/support/install/2021/Maple/Install.html


#### To run `Maple` benchmarks

1. Execute the following command in your favorite terminal from this directory

```
> maple maplescripts.mpl
```