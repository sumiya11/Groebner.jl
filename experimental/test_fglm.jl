if !isdefined(Main, :Groebner)
    using Groebner
end
using AbstractAlgebra

sys = Groebner.noonn(8, k=GF(2^26 + 15))

R, (x, y) = polynomial_ring(GF(2^26 + 15), ["x", "y"])
sys = [x * y + 1, y^200 + 1]

# gb_lex = Groebner.groebner(sys, ordering=Groebner.Lex())
gb_drl = Groebner.groebner(sys, ordering=Groebner.DegRevLex());

minpoly_last = Groebner.fglm_residuals_in_batch(
    gb_drl,
    Groebner.DegRevLex(),
    Groebner.Lex(),
    statistics=:timings,
    loglevel=0
)

gb_fglm = Groebner.fglm(gb_drl, Groebner.DegRevLex(), Groebner.Lex(), statistics=:timings);
