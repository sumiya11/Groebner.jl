module Groebner

debug() = false

"""
    always_vectorize()

Always try to vectorize hot loops using SIMD intrinsics.
"""
always_vectorize() = false

# For compatibility with polynomial types from AbstractAlgebra
import AbstractAlgebra
import AbstractAlgebra: base_ring, elem_type

import Combinatorics

import Logging
import Logging: ConsoleLogger, LogLevel

# For compatibility with polynomial types from MultivariatePolynomials
import MultivariatePolynomials
import MultivariatePolynomials: AbstractPolynomial, AbstractPolynomialLike

import Primes
import Primes: nextprime

import Random

# Some simple reference implementations
include("reference.jl")

# CRT and rational reconstruction
include("arithmetic/modular.jl")

# Supported monomial orderings
include("monoms/orderings.jl")
# Monomial implementations
include("monoms/packedutils.jl")
include("monoms/powervector.jl")
include("monoms/packedpairs.jl")

# Type aliases for internal objects
include("types.jl")

# Control of keywords in groebner and friends,
# and selection of optimal hyper-parameters
include("keywords.jl")

# Input-output conversions for polynomials
include("io.jl")

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
# All together combined in a generic f4 implementation 
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

export groebner
export isgroebner
export normalform

export fglm
export kbase

export Lex, DegLex, DegRevLex, InputOrdering
export NotPacked, Packed, best

# Set the contents of README.md as the docstring to this module
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Groebner

# 라헬
end
