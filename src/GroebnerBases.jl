
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
                        content, change_base_ring

import Primes: nextprime

# hmm??
import Nemo

using Combinatorics

using UnicodePlots


include("modular.jl")

include("common.jl")

include("testgens.jl")

include("buchberger.jl")

# f4 implementation over finite fields
include("f4.jl")

# fglm implementation over finite fields
include("fglm.jl")

# everything above is composed here
include("algorithm.jl")


export f4
export rootn
export groebner

end
