using DynamicPolynomials
using Groebner
import Profile

macro pr(ex, args...)
    return quote
        Profile.clear()
        Profile.init(n=10^7, delay=1e-5)
        Profile.@profile $(esc(ex))
        view_profile(; $(esc.(args)...))
    end
end

include((@__DIR__)*"/runge-kutta-6-6.jl");

vs = split("a_21 a_31 a_32 a_41 a_42 a_43 a_51 a_52 a_53 a_54 a_61 a_62 a_63 a_64 a_65 b_1 b_2 b_3 b_4 b_5 b_6")
length(vs)
length(system)

@time Groebner.groebner(system,linalg=:prob);
