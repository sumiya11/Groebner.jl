using StructuralIdentifiability
using AbstractAlgebra, BenchmarkTools, JLD2
import Nemo, Profile

macro myprof(ex)
    :((VSCodeServer.Profile).clear();
    Profile.init(n=10^7, delay=0.0000001);
    Profile.@profile $ex;
    VSCodeServer.view_profile(;))
end

ode = @ODEmodel(
    EGFR'(t) =
        EGFR_turnover * pro_EGFR(t) + EGF_EGFR(t) * reaction_1_k2 -
        EGFR(t) * EGFR_turnover - EGF_EGFR(t) * reaction_1_k1,
    pEGFR'(t) =
        EGF_EGFR(t) * reaction_9_k1 - pEGFR(t) * reaction_4_k1 +
        pEGFR_Akt(t) * reaction_2_k2 +
        pEGFR_Akt(t) * reaction_3_k1 - Akt(t) * pEGFR(t) * reaction_2_k1,
    pEGFR_Akt'(t) =
        Akt(t) * pEGFR(t) * reaction_2_k1 - pEGFR_Akt(t) * reaction_3_k1 -
        pEGFR_Akt(t) * reaction_2_k2,
    Akt'(t) =
        pAkt(t) * reaction_7_k1 + pEGFR_Akt(t) * reaction_2_k2 -
        Akt(t) * pEGFR(t) * reaction_2_k1,
    pAkt'(t) =
        pAkt_S6(t) * reaction_5_k2 - pAkt(t) * reaction_7_k1 +
        pAkt_S6(t) * reaction_6_k1 +
        pEGFR_Akt(t) * reaction_3_k1 - S6(t) * pAkt(t) * reaction_5_k1,
    S6'(t) =
        pAkt_S6(t) * reaction_5_k2 + pS6(t) * reaction_8_k1 -
        S6(t) * pAkt(t) * reaction_5_k1,
    pAkt_S6'(t) =
        S6(t) * pAkt(t) * reaction_5_k1 - pAkt_S6(t) * reaction_6_k1 -
        pAkt_S6(t) * reaction_5_k2,
    pS6'(t) = pAkt_S6(t) * reaction_6_k1 - pS6(t) * reaction_8_k1,
    EGF_EGFR'(t) =
        EGF_EGFR(t) * reaction_1_k1 - EGF_EGFR(t) * reaction_9_k1 -
        EGF_EGFR(t) * reaction_1_k2,
    y1(t) = a1 * (pEGFR(t) + pEGFR_Akt(t)),
    y2(t) = a2 * (pAkt(t) + pAkt_S6(t)),
    y3(t) = a3 * pS6(t)
)

siwr = @ODEmodel(
    S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
    W'(t) = xi * (I(t) - W(t)),
    R'(t) = gam * I(t) - (mu + a) * R(t),
    y(t) = k * I(t)
)

G = StructuralIdentifiability.ideal_generators(ode);

par = parent(G[1])
FF = Nemo.GF(2^31 - 1)
point = Nemo.QQ.(rand(1:100, length(AbstractAlgebra.gens(base_ring(base_ring(par))))))
G_zz = map(
    f -> map_coefficients(
        c -> evaluate(numerator(c), point) // evaluate(denominator(c), point),
        f
    ),
    G
);
G_zp = map(
    f -> map_coefficients(c -> FF(BigInt(numerator(c))) // FF(BigInt(denominator(c))), f),
    G_zz
);
typeof(parent(G_zp[1]))
typeof(G_zp[1])

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

gb = Groebner.groebner(G_zp);
@time graph, gb_1 = Groebner.groebner_learn(G_zp);
graph
graph.matrix_infos
flag, gb_2 = Groebner.groebner_apply!(graph, G_zp)
@assert flag && (gb == gb_1 == gb_2)

@btime Groebner.groebner($G_zp);
@btime Groebner.groebner_apply!($graph, $G_zp);

@myprof begin
    for i in 1:100
        Groebner.groebner_apply!(graph, G_zp)
    end
end

R, x = polynomial_ring(Nemo.GF(2^31 - 1), ["x$i" for i in 1:15], ordering=:degrevlex)

f = [a^rand(1:3) * b^rand(1:3) + c^rand(1:2) for a in x for b in x for c in x];

graph, gb = Groebner.groebner_learn(f);

@benchmark flag, gb_2 = Groebner.groebner_apply!($graph, $f)
@benchmark Groebner.groebner($f)

@myprof begin
    for i in 1:100
        Groebner.groebner_apply!(graph, f)
    end
end

c = Groebner.rootn(9, ground=Nemo.GF(2^31 - 1), ordering=:degrevlex)
graph, gb_1 = Groebner.groebner_learn(c);
flag, gb_2 = Groebner.groebner_apply!(graph, c);

@myprof begin
    for i in 1:100
        Groebner.groebner_apply!(graph, c)
    end
end

c = Groebner.rootn(9, ground=Nemo.GF(2^31 - 1), ordering=:degrevlex)
