
module Groebner

debug() = false

PRIMES = 0
F4TIME = 0.0
RECTIME = 0.0
CORRTIME = 0.0

ADD_ROW_HIT  = 0
ADD_ROW_MISS = 0
ADD_ROW_COLLISION = 0

function printhash()
    global ADD_ROW_HIT
    global ADD_ROW_MISS
    global ADD_ROW_COLLISION

    println("ADD_ROW_HIT = $(ADD_ROW_HIT)")
    println("ADD_ROW_MISS = $(ADD_ROW_MISS)")
    println("ADD_ROW_COLLISION = $(ADD_ROW_COLLISION)")

    ADD_ROW_HIT  = 0
    ADD_ROW_MISS = 0
    ADD_ROW_COLLISION = 0
end
ADD_ROW_HIT  = 0
ADD_ROW_MISS = 0

function printall()
    global PRIMES
    global F4TIME
    global RECTIME
    global CORRTIME
    s = sum((F4TIME, RECTIME, CORRTIME))
    println("run success,")
    println("PRIMES = $PRIMES")
    println("F4TIME = $(F4TIME/s)")
    println("RECTIME = $(RECTIME/s)")
    println("CORRTIME = $(CORRTIME/s)")
    PRIMES = 0
    F4TIME = 0.0
    RECTIME = 0.0
    CORRTIME = 0.0
end

# most of these are not needed
# TODO
import AbstractAlgebra
import AbstractAlgebra.Generic: MPoly, GFElem, MPolyRing, Poly
import AbstractAlgebra: leading_term, QQ, PolynomialRing, terms,
                        coeff, divides, base_ring, elem_type,
                        rref, isconstant, leading_coefficient,
                        map_coefficients, monomials, degree,
                        degrees, isconstant, leading_monomial,
                        GF, gens, MatrixSpace, coefficients,
                        crt, ordering, exponent_vectors, lift,
                        MPolyBuildCtx, finish, push_term!, ZZ,
                        content, change_base_ring, exponent_vector,
                        lcm, monomial, RingElem, set_exponent_vector!,
                        nvars, data, characteristic, isdivisible_by,
                        divexact, symbols

# for example systems
import Combinatorics

import Primes
import Primes: nextprime

import Random

import Logging
import Logging: ConsoleLogger, LogLevel

import MultivariatePolynomials
import MultivariatePolynomials: AbstractPolynomial, AbstractPolynomialLike

# type aliases for internal objects
include("internaltypes.jl")

# some simple reference implementations
include("reference.jl")

# computation parameters control
include("metainfo.jl")

# input-output conversions for polynomials
include("io.jl")

# modular arithmetic for CRT and rational reconstruction
include("arithmetic/modular.jl")

#= generic f4 implementation =#
#= the heart of this library =#
# `MonomialHashtable` implementation
include("f4/hashtable.jl")
# `Pairset` and `Basis` implementations
include("f4/basis.jl")
# `MacaulayMatrix` implementation
include("f4/matrix.jl")
include("f4/sorting.jl")
# Lucky prime numbers implementation
include("f4/lucky.jl")
# Buffers for manipulations with polynomial coefficients
include("f4/coeffs.jl")
# Enable tracing
include("f4/tracer.jl")
# All together combined in generic f4 implementation 
include("f4/f4.jl")

#= high level functions =#
include("gb/groebner.jl")
include("gb/isgroebner.jl")
include("gb/normalform.jl")
include("gb/correctness.jl")
# include("f4/statistics.jl")

#= generic fglm implementation =#
include("fglm/linear.jl")
include("fglm/fglm.jl")

# api
include("interface.jl")

# example systems definitions
include("testgens.jl")


export groebner
export isgroebner
export normalform

export fglm
export kbase

end
