module Groebner
# Groebner is a package for computing Gröbner bases.
#
# Groebner works over integers modulo a prime and over the rationals. Groebner
# does not provide a polynomial implementation of its own but relies on existing
# symbolic computation packages in Julia for communicating with the user
# instead. At its heart Groebner implements F4 and modular techniques.

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

using ExprTools: splitdef, combinedef

import Logging
import Logging: ConsoleLogger, LogLevel

import MultivariatePolynomials
import MultivariatePolynomials: AbstractPolynomial, AbstractPolynomialLike

import Primes
import Primes: nextprime

import Random

import Base.Threads: atomic_xchg!

include("utils/logging.jl")
include("utils/invariants.jl")

# Fast arithmetic in Z_p
include("arithmetic/Z_p.jl")
# Chinese remainder theorem and rational number reconstruction
include("arithmetic/reconstruction.jl")

# Supported monomial orderings
include("monoms/orderings.jl")
# Monomial implementations
include("monoms/packedutils.jl")
include("monoms/powervector.jl")
include("monoms/packedpairs.jl")

# type aliases for internal objects
# TODO(Alex): eliminate this file!
include("types.jl")

# control of keywords in groebner
# and optimal parameter selection
include("keywords.jl")

# input-output conversions for polynomials
include("input-output/common.jl")
include("input-output/AbstractAlgebra.jl")
include("input-output/DynamicPolynomials.jl")

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
