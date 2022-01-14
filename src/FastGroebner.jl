
module FastGroebner

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

# for testing
import Combinatorics

# controversial
using LoopVectorization

import Primes
import Primes: nextprime

import Random

# we want RadixSort to sort exponents (?)
import SortingAlgorithms

using PrettyPrinting

import Logging
import Logging: ConsoleLogger, LogLevel

#
include("utils.jl")

#
include("common.jl")

# input-output conversions for polynomials
include("io.jl")

#= operations with numbers =#
# unsigned arithmetic for finite field computations
include("arithmetic/unsigned.jl")
# modular arithmetic for CRT and rational reconstruction
include("arithmetic/modular.jl")
# basis coefficients manipulations
include("arithmetic/coeffs.jl")

#= f4 implementation over finite fields =#
include("f4/symbolic.jl")
include("f4/hash.jl")
include("f4/linear.jl")
include("f4/sorting.jl")
include("f4/f4.jl")
include("f4/normalform.jl")
# include("f4/statistics.jl")

#= fglm implementation over finite fields =#
include("fglm/fglm.jl")

# the heart of this library
include("groebner.jl")

# example systems definitions
include("testgens.jl")


export groebner
export isgroebner

end
