using DynamicPolynomials, Nemo
using Groebner
using Printf

include((@__DIR__) * "/rss_tracker.jl")
include((@__DIR__) * "/aa-runge-kutta-8-7_gf.jl");

setup_memuse_tracker()

const maxpairs = 500

# Compile
n = Groebner.noonn(3);
Groebner.groebner(n, ordering=DegRevLex(), maxpairs=maxpairs);

# Run!
@time gb = Groebner.groebner(
    system,
    ordering=DegRevLex(),
    maxpairs=maxpairs,
    monoms=Packed{UInt8}()
);

GC.gc()
