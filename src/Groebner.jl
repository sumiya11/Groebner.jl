module Groebner
# Groebner is a package for computing Gröbner bases. This is the main file.

# Groebner works over integers modulo a prime and over the rationals. At its
# heart Groebner implements F4 and modular techniques.

"""
    invariants_enabled() -> Bool

Specifies if custom asserts and invariants are checked.
If `false`, then all checks are disabled, and entail no runtime overhead.

See also `@invariant` in `src/utils/invariants.jl`.
"""
invariants_enabled() = true

"""
    logging_enabled() -> Bool

Specifies if custom logging is enabled.
If `false`, then all logging is disabled, and entails no runtime overhead.

See also `@log` in `src/utils/logging.jl`.
"""
logging_enabled() = true

# Groebner does not provide
# a polynomial implementation of its own but relies on existing symbolic
# computation packages in Julia for communicating with the user instead.
# Groebner accepts as its input polynomials from the Julia packages
# AbstractAlgebra.jl (Oscar.jl) and MultivariatePolynomials.jl. This list can be
# extended; see `src/input-output.jl` for details. TODO: make the
# AbstractAlgebra and MultivariatePolynomials dependencies optional!
import AbstractAlgebra
import AbstractAlgebra: base_ring, elem_type

import Combinatorics

using Logging

import MultivariatePolynomials
import MultivariatePolynomials: AbstractPolynomial, AbstractPolynomialLike

import Primes
import Primes: nextprime

import Random

import Base.Threads: Atomic, threadid, atomic_xchg!
import Base.MultiplicativeInverses: UnsignedMultiplicativeInverse

include("utils/logging.jl")
include("utils/invariants.jl")
include("utils/keywords.jl")
include("utils/testsystems.jl")

# Supported monomial orderings
include("monomials/orderings.jl")
# Monomial implementations
include("monomials/packedutils.jl")
include("monomials/exponentvector.jl")
include("monomials/packedvector.jl")

# Defines some type aliases for Groebner
include("utils/types.jl")

# Fast arithmetic in Z_p
include("arithmetic/Zp.jl")
# Not so fast arithmetic in the rationals
include("arithmetic/QQ.jl")

# Selecting algorithm parameters
include("utils/parameters.jl")

# Input-output conversions for polynomials
include("input-output/input-output.jl")
include("input-output/AbstractAlgebra.jl")
# include("input-output/DynamicPolynomials.jl")

#= generic f4 implementation =#
#= the heart of this library =#
# `MonomialHashtable` implementation
include("f4/hashtable.jl")
# `Pairset` and `Basis` implementations
include("f4/basis.jl")
# Computation graph for F4
include("f4/graph.jl")
# `MacaulayMatrix` implementation
include("f4/matrix.jl")
include("f4/sorting.jl")
# Some tracing
include("f4/tracer.jl")
# All together combined in F4 algorithm
include("f4/f4.jl")
include("f4/learn-apply.jl")

#= more high level functions =#
# Lucky prime numbers
include("groebner/lucky.jl")
# Rational number reconstruction and CRT reconstruction
include("groebner/reconstruction.jl")
# GB state
include("groebner/state.jl")
# Correctness checks
include("groebner/correctness.jl")
# `groebner` backend
include("groebner/groebner.jl")
# `isgroebner` backend
include("groebner/isgroebner.jl")
# `normalform` backend
include("groebner/normalform.jl")

#= generic fglm implementation =#
include("fglm/linear.jl")
include("fglm/fglm.jl")

# API
include("interface.jl")

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
