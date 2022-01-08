
module GroebnerBases

import AbstractAlgebra
import AbstractAlgebra.Generic: MPoly, GFElem
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
                        MPolyRing, nvars, data, characteristic

import Combinatorics
import Primes
import Primes: nextprime
import SparseArrays: SparseVector, findnz, nnz
import Random

import DataStructures: SortedSet

# we want RadixSort to sort exponents (?)
import SortingAlgorithms

import IterTools

using PrettyPrinting

import Logging
import Logging: ConsoleLogger, LogLevel

import Cthulhu

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
include("f4/statistics.jl")


# the heart of this library
include("groebner.jl")

# example systems definitions
include("testgens.jl")

export groebner
export isgroebner

end
