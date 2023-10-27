@def title = "Groebner.jl — For developers"
@def hasmath = true
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->

## Development

```julia:load_groebner
using Groebner # hide
```

### Project source

You can get a clone of the GitHub repository with:

```
git clone https://github.com/sumiya11/Groebner.jl
```

If you want to browse the repository online, or fork it on GitHub, it can be accessed here:

[https://github.com/sumiya11/Groebner.jl](https://github.com/sumiya11/Groebner.jl)

The up-to-date branch in the repository is "master".

### Development directions

WIP!

### Utilities in Groebner.jl

Internally, Groebner.jl provides utilities for convenient assertion checking, logging, and timings measurment. These are implemented via macros `@invariant`, `@log`, and `@timeit`.

Each utility can be globally enabled/disabled by setting the appropriate switch to `true`/`false`. The switches and their default values are: 

| Macro        | Switch                         | Default value |
|--------------|--------------------------------|---------------|
| `@invariant` | `invariants_enabled`           | `true`        |
| `@log`       | `logging_enabled`              | `true`        |
| `@timeit`    | `performance_counters_enabled` | `false`       |

If the switch value is `false`, then the corresponding macro gets *compiled-out*, and has no performance impact.

Therefore, for example, if you would like to get the best performance out of Groebner.jl, and do not care about assertions and logging, you can do

```julia
using Groebner
Groebner.invariants_enabled() = false
Groebner.logging_enabled() = false
```

#### Logging

Assuming `Groebner.logging_enabled()` is `true`, for printing some logs you can do, for example,
```julia:./dev-log
using Groebner, AbstractAlgebra

R, (x,y,z) = GF(2^31-1)["x","y","z"]

# Will print some debug info.
# Use lower loglevel, e.g., loglevel=-5, to print A LOT of info
groebner([x*y + z, x*z + y], loglevel=-2)
```

#### Measuring performance

All functions in the interface have keyword argument `statistics`. This argument can be set to either of: `:no`, `:timings`, and `:all`. 

Use `statistics=:timings` to print the table with timings and allocations of internal functions of Groebner.jl.

Since `Groebner.performance_counters_enabled()` is `false` by default, you should set it to `true` to record runtime statistics. For example,

```julia:dev-timings
using Groebner, AbstractAlgebra

# Enable performance counters
Groebner.performance_counters_enabled() = true

R, (x,y,z) = GF(2^31-1)["x","y","z"]

groebner([x*y + z, x*z + y], statistics=:no) #hide
groebner([x*y + z, x*z + y], statistics=:timings)
```
