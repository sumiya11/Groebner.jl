module Groebner
# Groebner is a package for computing Gröbner bases.

# Groebner works over integers modulo a prime and over the rationals. At its
# heart Groebner implements F4 and modular techniques. Groebner does not provide
# a polynomial implementation of its own but relies on existing symbolic
# computation packages in Julia for communicating with the user instead.

"""
    invariants_enabled()

Specifies if custom asserts and invariants are checked.
If `false`, then all checks are disabled, and entail no runtime overhead.

See also `enable_invariants`.
"""
invariants_enabled() = true

"""
    logging_enabled()

Specifies if custom logging is enabled.
If `false`, then all logging is disabled, and entails no runtime overhead.

See also `enable_logging`.
"""
logging_enabled() = true

# Groebner accepts as input polynomials from the Julia packages
# AbstractAlgebra.jl (Oscar.jl) and MultivariatePolynomials.jl.
# This list can be extended; see `src/input-output.jl` for details.
# TODO: make the AbstractAlgebra and MultivariatePolynomials dependencies optional!
import AbstractAlgebra
import AbstractAlgebra: base_ring, elem_type

import Combinatorics

import MultivariatePolynomials
import MultivariatePolynomials: AbstractPolynomial, AbstractPolynomialLike

import Primes
import Primes: nextprime

import Random

import Base.Threads: Atomic, threadid, atomic_xchg!
 
include("utils/logging.jl")
include("utils/invariants.jl")
include("utils/keywords.jl")

# Supported monomial orderings
include("monomials/orderings.jl")
# Monomial implementations
include("monomials/packedutils.jl")
include("monomials/powervector.jl")
include("monomials/packedpairs.jl")

# Fast arithmetic in Z_p
include("arithmetic/Z_p.jl")

# Defines some type aliases for Groebner
include("utils/types.jl")

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
# `MacaulayMatrix` implementation
include("f4/matrix.jl")
include("f4/sorting.jl")
# Enable tracing
include("f4/tracer.jl")
# All together combined in F4 algorithm
include("f4/f4.jl")

#= more high level functions =#
# Lucky prime numbers
include("groebner/lucky.jl")
# Manipulation with big polynomial coefficients
include("groebner/coeffs.jl")
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
