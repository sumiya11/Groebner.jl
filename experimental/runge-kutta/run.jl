using DynamicPolynomials
using Groebner
import Profile
import Nemo

macro pr(ex, args...)
    return quote
        Profile.clear()
        Profile.init(n=10^8, delay=1e-5)
        Profile.@profile $(esc(ex))
        view_profile(; $(esc.(args)...))
    end
end

using Nemo
include((@__DIR__) * "/aa-runge-kutta-8-7.jl");

@time gb = Groebner.groebner(system, ordering=Groebner.DegRevLex(), maxpairs=1000);

vs = split(
    "a_21 a_31 a_32 a_41 a_42 a_43 a_51 a_52 a_53 a_54 a_61 a_62 a_63 a_64 a_65 a_71 a_72 a_73 a_74 a_75 a_76 a_81 a_82 a_83 a_84 a_85 a_86 a_87 b_1 b_2 b_3 b_4 b_5 b_6 b_7 b_8"
)
vs
ord = Groebner.WeightedOrdering([[10 for i in 1:28]..., [1 for i in 1:8]...])

gb = Groebner.groebner(
    system,
    linalg=:prob
    # ordering=ord
);

println(gb)

using AbstractAlgebra
include((@__DIR__) * "/aa-runge-kutta-6-6.jl");
length(system)
gens(parent(system[1]))

system_autoreduced = deepcopy(system);
for i in 1:length(system)
    f = system_autoreduced[i]
    I = system[i .!= 1:length(system)]
    system_autoreduced[i] = AbstractAlgebra.normal_form(f, I)
end;

@pr gb = Groebner.groebner(system, linalg=:prob, ordering=DegRevLex());

@polyvar x y z a b c
s = [y, x, z, a + b^2 + c^3]
g = groebner(s, ordering=Groebner.WeightedOrdering([1, 10, 10, 1, 1, 1]))
