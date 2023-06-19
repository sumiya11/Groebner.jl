module Groebner
# Groebner is a package for computing Gröbner bases.
#
# Groebner works over integers modulo a prime and over the rationals. Groebner
# does not provide a polynomial implementation of its own but relies on existing
# symbolic computation packages in Julia for communicating with the user
# instead. At its heart Groebner implements F4 and modular techniques.
#
# Here is how Groebner works from bird's eye view. 
# The primary function exported by Groebner is `groebner`.
#
#     (1)             (2)             (3)            (4)            (5)
#   polynoms  -->  internal  -->  F4 backend  -->  internal  -->  polynoms
#                    repr.                           repr.
#
# (1) -> (2): The `groebner` function accepts an array of polynomials and
# extracts the raw data: polynomial ring information, polynomial exponent
# vectors, polynomial coefficients. Then, it converts the data to an internal
# representation, which may depend on the features of input. For example,
# exponent vectors could be stored densely or sparsely. See
# `src/input-output.jl` for details. 
# (2) -> (3): The appropriate algorithm is called (Z_p vs. Q). The algorithm
# initializes the data structures necessary for F4 (hashtables for monomials,
# storages for critical pairs and polynomials) and invokes the F4 backend.
# Over the rationals modular computations are used to reduce the task to Z_p.
# This logic is implemented in `src/gb`. The F4 backend is encapsulated within
# the `src/f4` directory.
# (3) -> (4), (4) -> (5): The correctness of result is checked and some
# post-processing steps may be applied. Finally, the internal representation is
# extracted from the F4 data structures and converted back to the original
# polynomial types.

# Adapted from 
#   https://discourse.julialang.org/t/assert-alternatives/24775
"""
    enable_invariants()

Enforce runtime invariant checking. If `true`, checks are enabled. If
`false`, checks are disabled and entail no runtime overhead.
"""
enable_invariants() = true

"""
    @invariant expr

Throw an `AssertionError` if `expr` evaluates to `false`.
Can be enabled/disabled by `enable_invariants`.

**Note:** the Julia's `@assert` will trigger independently of the value returned by `enable_invariants`.
"""
macro invariant(expr)
    esc(:(
        if $(@__MODULE__).enable_invariants()
            @assert $expr
        end
    ))
end

macro record(condition, expr)
    esc(:(

    ))
end

# Groebner accepts as an input polynomials from the Julia packages
# AbstractAlgebra.jl (Oscar.jl) and MultivariatePolynomials.jl.
# This list can be extended; see `src/input-output.jl` for details.
import AbstractAlgebra
import AbstractAlgebra: base_ring, elem_type

import Combinatorics

import Logging
import Logging: ConsoleLogger, LogLevel

import MultivariatePolynomials
import MultivariatePolynomials: AbstractPolynomial, AbstractPolynomialLike

import Primes
import Primes: nextprime

import Random

# CRT and rational number reconstruction
include("arithmetic/reconstruction.jl")
# arithmetic in Z_p
include("arithmetic/modular.jl")

# supported monomial orderings
include("monoms/orderings.jl")
# monomial implementations
include("monoms/packedutils.jl")
include("monoms/powervector.jl")
include("monoms/packedpairs.jl")

# type aliases for internal objects
include("types.jl")

# control of keywords in groebner
# and optimal parameter selection
include("keywords.jl")

# input-output conversions for polynomials
include("input-output.jl")

#= generic f4 implementation =#
#= the heart of this library =#
# `MonomialHashtable` implementation
include("f4/hashtable.jl")
# `Pairset` and `Basis` implementations
include("f4/basis.jl")
# `MacaulayMatrix` implementation
include("f4/matrix.jl")
include("f4/sorting.jl")
# Enable tracing
include("f4/tracer.jl")
# All together combined in F4 algorithm
include("f4/f4.jl")

#= more high level functions =#
# Lucky prime numbers
include("gb/lucky.jl")
# Manipulation with big polynomial coefficients
include("gb/coeffs.jl")
# Correctness checks
include("gb/correctness.jl")
# `groebner` backend
include("gb/groebner.jl")
# `isgroebner` backend
include("gb/isgroebner.jl")
# `normalform` backend
include("gb/normalform.jl")

#= generic fglm implementation =#
include("fglm/linear.jl")
# `fglm` and `kbase` backend
include("fglm/fglm.jl")

# api
include("interface.jl")

# example systems definitions
include("testsystems.jl")

using SnoopPrecompile
include("precompile.jl")

export groebner
export isgroebner
export normalform

export kbase

export Lex,
    DegLex, DegRevLex, InputOrdering, WeightedOrdering, BlockOrdering, MatrixOrdering
export NotPacked, Packed, best_monom_representation

@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Groebner

# 라헬
end # module Groebner
