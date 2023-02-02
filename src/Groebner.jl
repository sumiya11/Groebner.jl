module Groebner

debug() = false

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

using Base.Threads

include("utils.jl")

# some simple reference implementations
include("reference.jl")

# CRT and rational reconstruction
include("arithmetic/modular.jl")
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
# All together combined in generic f4 implementation 
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

# 라헬
end
