using Nemo
# using Groebner
using Printf

include((@__DIR__) * "/rss_tracker.jl")
include((@__DIR__) * "/aa-runge-kutta-8-7_gf.jl");

setup_memuse_tracker()

# Run!
@time gb = Groebner.groebner(system, ordering=Groebner.DegRevLex(), loglevel=:debug);
