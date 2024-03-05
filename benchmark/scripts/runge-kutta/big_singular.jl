import Pkg
Pkg.activate("../singular")

import AbstractAlgebra
import Singular
using Nemo

include((@__DIR__) * "/rss_tracker.jl")
include((@__DIR__) * "/aa-runge-kutta-6-6.jl");

function aa_to_singular(poly)
    Rxx = parent(poly)
    Rqq = AbstractAlgebra.base_ring(Rxx)
    xstrings = map(string, AbstractAlgebra.gens(Rxx))
    base = Singular.QQ
    new_ring, _ = Singular.polynomial_ring(base, xstrings, internal_ordering=:degrevlex)
    AbstractAlgebra.change_base_ring(
        AbstractAlgebra.base_ring(new_ring),
        poly,
        parent=new_ring
    )
end

setup_memuse_tracker()

singular_system = map(aa_to_singular, system)
singular_ring = parent(singular_system[1])

println(ordering(singular_ring))

# Run!
singular_ideal = Singular.Ideal(singular_ring, singular_system)
# @time gb = Singular.std(singular_ideal);

# Run slimgb
singular_ideal = Singular.Ideal(singular_ring, singular_system)
@time gb = Singular.slimgb(singular_ideal);

GC.gc()
