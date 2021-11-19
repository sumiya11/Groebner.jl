
module GroebnerBases

import AbstractAlgebra
import AbstractAlgebra.Generic: MPoly, GFElem
import AbstractAlgebra: leading_term, QQ, PolynomialRing, terms,
                        coeff, divides, base_ring, elem_type,
                        rref, isconstant, leading_coefficient,
                        map_coefficients, monomials, degree,
                        degrees, isconstant, leading_monomial,
                        GF, gens, MatrixSpace

using Combinatorics

using UnicodePlots


include("modular.jl")

include("common.jl")

include("testgens.jl")

include("buchberger.jl")

# main functionality
include("f4.jl")


export f4
export rootn

end
