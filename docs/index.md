@def title = "SymbolicUtils.jl — Symbolic programming in Julia"
@def hasmath = false
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->


# Groebner


`Groebner.jl` is a package for fast and generic Gröbner bases computations
based on Faugère's F4 algorithm [[1]](https://www-polsys.lip6.fr/~jcf/Papers/F99a.pdf).

## Installation

To install Groebner.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("https://github.com/sumiya11/Groebner.jl")
```

## Interface

```julia:load_groebner
using Groebner # hide
```

{{doc groebner groebner fn}}

{{doc isgroebner isgroebner fn}}

{{doc reducegb reducegb fn}}

