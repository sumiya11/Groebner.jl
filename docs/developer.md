@def title = "Groebner.jl — For developers"
@def hasmath = true
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->

## For developers

```julia:load_groebner
using Groebner # hide
```

### Getting Groebner.jl for development

git clone ...

or Pkg.dev

### Structure of Groebner.jl

..Graph of function calls..

### Utilities

#### Assertions

To enable/disable all `@invariant` macro in Groebner.jl, change `Groebner.invariants_enabled()` to `true`/`false`. 

{{doc Groebner.@invariant Groebner.@invariant macro}}

#### Logging

To enable/disable all `@log` macro in Groebner.jl, change `Groebner.logging_enabled()` to `true`/`false`. 

For example,
```julia:dev-log
using Groebner, AbstractAlgebra
_, (x,y,z) = GF(2^31-1)["x","y","z"]
groebner([x*y + z, x*z + y], loglevel=-3)
```

{{doc Groebner.@log Groebner.@log macro}}

#### Collecting statistics and measuring performance

All functions in the interface of Groebner.jl have keyword argument
`statistics`. This keyword argument can be set to either of: `:no`, `:timings`, and `:all`. 

To enable/disable all `@log` macro in Groebner.jl, change `Groebner.logging_enabled()` to `true`/`false`. 

For example,
```julia:dev-timings
using Groebner, AbstractAlgebra
_, (x,y,z) = GF(2^31-1)["x","y","z"]
groebner([x*y + z, x*z + y], statistics=:timings)
```

Timings are recorded

{{doc Groebner.@timeit Groebner.@timeit macro}}
