
module Groebner

# HEURISTIC = 0
# RANDOMIZED = 0

# most of these are not needed
# TODO
import AbstractAlgebra
import AbstractAlgebra.Generic: MPoly, GFElem, MPolyRing
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
                        divexact

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

# so, what is that exactly?
include("utils.jl")

# some simple reference implementations
include("common.jl")

# input-output conversions for polynomials
include("io.jl")

#= operations with numbers =#
# unsigned arithmetic for finite field computations
include("arithmetic/unsigned.jl")
# modular arithmetic for CRT and rational reconstruction
include("arithmetic/modular.jl")

#= f4 implementation over finite fields =#
# the heart of this library
include("f4/structs.jl")
include("f4/symbolic.jl")
include("f4/hash.jl")
include("f4/linear.jl")
include("f4/sorting.jl")
include("f4/lucky.jl")
include("f4/coeffs.jl")
include("f4/f4.jl")
include("f4/groebner.jl")
include("f4/isgroebner.jl")
include("f4/normalform.jl")
include("f4/correctness.jl")
# include("f4/statistics.jl")

#= fglm implementation over finite fields =#
include("fglm/fglm.jl")

# api
include("interface.jl")

# example systems definitions
include("testgens.jl")


export groebner
export isgroebner
export normalform

end
