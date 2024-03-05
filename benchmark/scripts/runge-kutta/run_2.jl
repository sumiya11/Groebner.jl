# using DynamicPolynomials, Nemo
# using Groebner

# include((@__DIR__) * "/rss_tracker.jl")
include((@__DIR__) * "/aa-runge-kutta-6-6.jl");

# setup_memuse_tracker()

Groebner.logging_enabled() = true

@time gb = Groebner.groebner(system, loglevel=-1, ordering=Groebner.DegRevLex());
