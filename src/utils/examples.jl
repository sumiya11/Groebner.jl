# This file is a part of Groebner.jl. License is GNU GPL v2.
module Examples
# Systems used for testing and benchmarking.
# Useful references:
#   - https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl
#   - https://gitlab.lip6.fr/eder/msolve-examples
#   - https://web.archive.org/web/20201202185136/http://www.cecm.sfu.ca/%7Erpearcea/mgb.html
#   - https://github.com/symbolicdata/data

using ..Groebner.AbstractAlgebra
using ..Groebner.Combinatorics
using ..Groebner.Random

#! format: off
# Syntax formatting is off in this file.

# Göran Björck and Ralf Fröberg: `A faster way to count the solutions of
# inhomogeneous systems of algebraic equations, with applications to cyclic
# n-roots', in J. Symbolic Computation (1991) 12, pp 329–336.
#
# source: https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/v0.1.6/src/systems.jl
#
# cyclic(8)
# Critical pair degree sequence: 2,3,4,5,6,7,8,9,10,11,12,5,6,7,8,9,10,11,12,13,6,7,8,9,10,11,12,13,6,7,8,9,10,11,12,13,14,15
function cyclicn(n; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, z = np.polynomial_ring(k, ["z$i" for i in 1:n], internal_ordering=internal_ordering)
    [
        (
            sum(prod(z[(k - 1) % n + 1] for k in j:(j + m)) for j in 1:n) for m in 0:(n - 2)
        )...,
        prod(z) - 1
    ]
end

# Not to be confused with cyclic-n !!
function rootn(n; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, xs = np.polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=internal_ordering)
    ans = [sum(map(prod, Combinatorics.combinations(xs, i))) for i in 1:n]
    ans[end] -= (-1)^(n - 1)
    ans
end

# reimern(7)
# Critical pair degree sequence: 3,4,5,6,7,8,9,8,7,7,6,6,5,6,5,5,5,4,4,3,4,5,6,7,8,9,10,11,12,13,14,15
function reimern(n; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, xs = np.polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=internal_ordering)
    [sum((-1)^(i + 1) * 2 * xs[i]^j for i in 1:n) - 1 for j in 2:(n + 1)]
end

# S. Katsura, W. Fukuda, S. Inawashiro, N.M. Fujiki and R. Gebauer,
#  Cell Biophysics, Vol 11, pages 309–319, 1987.
#
# source: https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/v0.1.6/src/systems.jl
#
# katsuran(10):
# Critical pair degree sequence: 2,3,4,5,6,7,8,9,10,11,12
function katsuran(n; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, x = np.polynomial_ring(k, ["x$i" for i in 0:n], internal_ordering=internal_ordering)
    [
        (
            sum(x[abs(l) + 1] * x[abs(m - l) + 1] for l = (-n):n if abs(m - l) <= n) -
            x[m + 1] for m in 0:(n - 1)
        )...,
        x[1] + 2sum(x[i + 1] for i in 1:n) - 1
    ]
end

# noonn(8):
# Critical pair degree sequence: 4,5,6,7,8,9,10,11,12,13,14,15,16,15,16
function noonn(n; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    without(x, k) = x[1:end .!= k]

    R, xs = np.polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=internal_ordering)
    fs = zeros(R, n)
    for i in 1:n
        other = without(xs, i)
        fs[i] = xs[i] * (10 * sum(other .^ 2) - 11) + 10
    end
    fs
end

###

# Critical pair degree sequence: 4,4,4,5,5,4,4,4,4,5,6
function hexapod(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (t1,t2,t3,a,b,c) = np.polynomial_ring(k, ["t1","t2","t3","a", "b", "c"], internal_ordering=internal_ordering)
    [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976] 
end

###
# Henrion

# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion5.ms
function henrion5(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (f1, f2, f3, f4, f5, t) =
        np.polynomial_ring(k, ["f1", "f2", "f3", "f4", "f5", "t"], internal_ordering=internal_ordering)
    [
        2*f1*f2*f3*f4*f5-9823275,
        k(21) // k(5)*f1*f2*f4*f5+k(16) // k(5)*f1*f3*f4*f5+k(9) // k(5)*f2*f3*f4*f5+k(24) // k(5)*f1*f2*f3*f5+5*f4*f3*f1*f2-4465125,
        k(14) // k(5)*f4*f5*f1+k(14) // k(5)*f4*f5*f2+k(8) // k(5)*f3*f4*f5+k(18) // k(5)*f1*f2*f5+k(24) // k(5)*f1*f3*f5+k(18) // k(5)*f2*f3*f5+4*f3*f1*f2+6*f1*f2*f4+6*f3*f4*f1+4*f2*f3*f4-441486,
        k(7) // k(5)*f4*f5+k(12) // k(5)*f5*f1+k(12) // k(5)*f5*f2+k(12) // k(5)*f5*f3+3*f1*f2+4*f3*f1+4*f4*f1+3*f2*f3+4*f4*f2+3*f3*f4-15498,
        k(6) // k(5)*f5+2*f4+2*f3+2*f2+2*f1-215,
        f1+2*f2+3*f3+4*f4+5*f5+6*t
    ]
end

# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion6.ms
function henrion6(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (f1, f2, f3, f4, f5, f6) =
        np.polynomial_ring(k, ["f1", "f2", "f3", "f4", "f5", "f6"], internal_ordering=internal_ordering)
    [
        2*f1*f2*f3*f4*f5*f6-1404728325,
        6*f5*f4*f3*f1*f2+k(11) // k(6)*f2*f3*f4*f5*f6+k(16) // k(3)*f1*f2*f3*f5*f6+k(9) // k(2)*f1*f2*f4*f5*f6+k(10) // k(3)*f1*f3*f4*f5*f6+k(35) // k(6)*f1*f2*f3*f4*f6-648336150,
        5*f4*f3*f1*f2+5*f2*f3*f4*f5+k(5) // k(3)*f3*f4*f5*f6+8*f1*f2*f3*f5+9*f1*f2*f4*f5+8*f1*f3*f4*f5+4*f1*f2*f5*f6+k(16) // k(3)*f1*f3*f5*f6+3*f1*f4*f5*f6+4*f2*f3*f5*f6+3*f2*f4*f5*f6+k(14) // k(3)*f1*f2*f3*f6+7*f1*f2*f4*f6+7*f1*f3*f4*f6+k(14) // k(3)*f2*f3*f4*f6-67597623,
        6*f1*f2*f5+8*f1*f3*f5+6*f2*f3*f5+k(8) // k(3)*f5*f6*f3+k(8) // k(3)*f5*f6*f2+k(8) // k(3)*f5*f6*f1+k(7) // k(2)*f1*f2*f6+k(14) // k(3)*f1*f3*f6+k(14) // k(3)*f1*f4*f6+k(7) // k(2)*f2*f3*f6+k(14) // k(3)*f2*f4*f6+k(7) // k(2)*f3*f4*f6+6*f4*f5*f1+k(3) // k(2)*f4*f5*f6+4*f3*f1*f2+4*f2*f3*f4+6*f3*f4*f1+4*f3*f4*f5+6*f1*f2*f4+6*f4*f5*f2-2657700,
        k(4) // k(3)*f5*f6+k(7) // k(3)*f6*f1+k(7) // k(3)*f6*f2+k(7) // k(3)*f6*f3+k(7) // k(3)*f6*f4+3*f1*f2+4*f3*f1+4*f4*f1+4*f5*f1+3*f2*f3+4*f4*f2+4*f5*f2+3*f3*f4+4*f5*f3+3*f4*f5-46243,
        k(7) // k(6)*f6+2*f5+2*f4+2*f3+2*f2+2*f1-358
    ]
end

# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion7.ms
#
# Critical pair degree sequence: 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
function henrion7(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (f1, f2, f3, f4, f5, f6, f7) = np.polynomial_ring(
        k,
        ["f1", "f2", "f3", "f4", "f5", "f6", "f7"],
        internal_ordering=internal_ordering
    )
    [
        2*f1*f2*f3*f4*f5*f6*f7-273922023375,
        k(45) // k(7)*f1*f2*f3*f4*f6*f7+k(40) // k(7)*f1*f2*f3*f5*f6*f7+k(33) // k(7)*f1*f2*f4*f5*f6*f7+k(24) // k(7)*f1*f3*f4*f5*f6*f7+k(48) // k(7)*f1*f2*f3*f4*f5*f7+7*f6*f5*f4*f3*f1*f2+k(13) // k(7)*f2*f3*f4*f5*f6*f7-127830277575,
        6*f5*f4*f3*f1*f2+6*f2*f3*f4*f5*f6+k(12) // k(7)*f3*f4*f5*f6*f7+10*f1*f2*f3*f4*f6+12*f1*f2*f3*f5*f6+12*f1*f2*f4*f5*f6+10*f1*f3*f4*f5*f6+k(36) // k(7)*f1*f2*f3*f6*f7+k(54) // k(7)*f1*f2*f4*f6*f7+k(30) // k(7)*f1*f2*f5*f6*f7+k(54) // k(7)*f1*f3*f4*f6*f7+k(40) // k(7)*f1*f3*f5*f6*f7+k(22) // k(7)*f1*f4*f5*f6*f7+k(36) // k(7)*f2*f3*f4*f6*f7+k(30) // k(7)*f2*f3*f5*f6*f7+k(22) // k(7)*f2*f4*f5*f6*f7+k(40) // k(7)*f1*f2*f3*f4*f7+k(64) // k(7)*f1*f2*f3*f5*f7+k(72) // k(7)*f1*f2*f4*f5*f7+k(64) // k(7)*f1*f3*f4*f5*f7+k(40) // k(7)*f2*f3*f4*f5*f7-13829872635,
        -585849123+5*f4*f3*f1*f2+5*f2*f3*f4*f5+5*f3*f4*f5*f6+8*f1*f2*f3*f5+9*f1*f2*f4*f5+8*f1*f3*f4*f5+k(11) // k(7)*f4*f5*f6*f7+8*f1*f2*f3*f6+12*f1*f2*f4*f6+9*f1*f2*f5*f6+12*f1*f3*f4*f6+12*f1*f3*f5*f6+8*f1*f4*f5*f6+8*f2*f3*f4*f6+9*f2*f3*f5*f6+8*f2*f4*f5*f6+k(27) // k(7)*f1*f2*f6*f7+k(36) // k(7)*f1*f3*f6*f7+k(36) // k(7)*f1*f4*f6*f7+k(20) // k(7)*f1*f5*f6*f7+k(27) // k(7)*f2*f3*f6*f7+k(36) // k(7)*f2*f4*f6*f7+k(20) // k(7)*f2*f5*f6*f7+k(27) // k(7)*f3*f4*f6*f7+k(20) // k(7)*f3*f5*f6*f7+k(32) // k(7)*f1*f2*f3*f7+k(48) // k(7)*f1*f2*f4*f7+k(48) // k(7)*f1*f2*f5*f7+k(48) // k(7)*f1*f3*f4*f7+k(64) // k(7)*f1*f3*f5*f7+k(48) // k(7)*f1*f4*f5*f7+k(32) // k(7)*f2*f3*f4*f7+k(48) // k(7)*f2*f3*f5*f7+k(48) // k(7)*f2*f4*f5*f7+k(32) // k(7)*f3*f4*f5*f7,
        -11675085+6*f1*f2*f6+8*f1*f3*f6+8*f1*f4*f6+6*f2*f3*f6+8*f2*f4*f6+6*f3*f4*f6+k(18) // k(7)*f6*f7*f4+k(18) // k(7)*f6*f7*f3+k(18) // k(7)*f6*f7*f2+k(18) // k(7)*f6*f7*f1+k(24) // k(7)*f1*f2*f7+k(32) // k(7)*f1*f3*f7+k(32) // k(7)*f1*f4*f7+k(32) // k(7)*f1*f5*f7+k(24) // k(7)*f2*f3*f7+k(32) // k(7)*f2*f4*f7+k(32) // k(7)*f2*f5*f7+k(24) // k(7)*f3*f4*f7+k(32) // k(7)*f3*f5*f7+k(24) // k(7)*f4*f5*f7+k(10) // k(7)*f5*f6*f7+6*f1*f2*f5+8*f1*f3*f5+6*f2*f3*f5+6*f5*f6*f3+6*f5*f6*f2+6*f5*f6*f1+6*f1*f2*f4+6*f4*f5*f2+6*f4*f5*f1+4*f4*f5*f6+6*f3*f4*f1+4*f3*f4*f5+4*f3*f1*f2+4*f2*f3*f4,
        k(9) // k(7)*f6*f7+k(16) // k(7)*f7*f1+k(16) // k(7)*f7*f2+k(16) // k(7)*f7*f3+k(16) // k(7)*f7*f4+k(16) // k(7)*f7*f5+3*f1*f2+4*f3*f1+4*f4*f1+4*f5*f1+4*f6*f1+3*f2*f3+4*f4*f2+4*f5*f2+4*f6*f2+3*f3*f4+4*f5*f3+4*f6*f3+3*f4*f5+4*f6*f4+3*f5*f6-116053,
        k(8) // k(7)*f7+2*f6+2*f5+2*f4+2*f3+2*f2+2*f1-553
    ]
end

# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion8.ms
function henrion8(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (f1, f2, f3, f4, f5, f6, f7, f8) = np.polynomial_ring(
        k,
        ["f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8"],
        internal_ordering=internal_ordering
    )
    [
        k(2)*f1*f2*f3*f4*f5*f6*f7*f8-k(69850115960625),
        k(8)*f7*f6*f5*f4*f3*f1*f2+k(15)//k(8)*f2*f3*f4*f5*f6*f7*f8+k(15)//k(2)*f1*f2*f3*f4*f5*f7*f8+k(55)//k(8)*f1*f2*f3*f4*f6*f7*f8+k(6)*f1*f2*f3*f5*f6*f7*f8+k(39)//k(8)*f1*f2*f4*f5*f6*f7*f8+k(7)//k(2)*f1*f3*f4*f5*f6*f7*f8+k(63)//k(8)*f1*f2*f3*f4*f5*f6*f8-k(32870642805000),
        -k(3654447799500)+k(27)//k(2)*f1*f2*f4*f5*f6*f8+k(45)//k(4)*f1*f3*f4*f5*f6*f8+k(27)//k(4)*f2*f3*f4*f5*f6*f8+k(7)//k(4)*f3*f4*f5*f6*f7*f8+k(12)*f1*f2*f3*f4*f5*f7+k(15)*f1*f2*f3*f4*f6*f7+k(16)*f1*f2*f3*f5*f6*f7+k(15)*f1*f2*f4*f5*f6*f7+k(12)*f1*f3*f4*f5*f6*f7+k(25)//k(4)*f1*f2*f3*f4*f7*f8+k(10)*f1*f2*f3*f5*f7*f8+k(11)//k(2)*f1*f2*f3*f6*f7*f8+k(45)//k(4)*f1*f2*f4*f5*f7*f8+k(33)//k(4)*f1*f2*f4*f6*f7*f8+k(9)//k(2)*f1*f2*f5*f6*f7*f8+k(10)*f1*f3*f4*f5*f7*f8+k(33)//k(4)*f1*f3*f4*f6*f7*f8+k(6)*f1*f3*f5*f6*f7*f8+k(13)//k(4)*f1*f4*f5*f6*f7*f8+k(25)//k(4)*f2*f3*f4*f5*f7*f8+k(11)//k(2)*f2*f3*f4*f6*f7*f8+k(9)//k(2)*f2*f3*f5*f6*f7*f8+k(13)//k(4)*f2*f4*f5*f6*f7*f8+k(27)//k(4)*f1*f2*f3*f4*f5*f8+k(45)//k(4)*f1*f2*f3*f4*f6*f8+k(27)//k(2)*f1*f2*f3*f5*f6*f8+k(7)*f6*f5*f4*f3*f1*f2+k(7)*f2*f3*f4*f5*f6*f7,
        -k(163221399000)+k(6)*f5*f4*f3*f1*f2+k(6)*f2*f3*f4*f5*f6+k(6)*f3*f4*f5*f6*f7+k(10)*f1*f2*f3*f4*f6+k(12)*f1*f2*f3*f5*f6+k(12)*f1*f2*f4*f5*f6+k(10)*f1*f3*f4*f5*f6+k(13)//k(8)*f4*f5*f6*f7*f8+k(10)*f1*f2*f3*f4*f7+k(16)*f1*f2*f3*f5*f7+k(12)*f1*f2*f3*f6*f7+k(18)*f1*f2*f4*f5*f7+k(18)*f1*f2*f4*f6*f7+k(12)*f1*f2*f5*f6*f7+k(16)*f1*f3*f4*f5*f7+k(18)*f1*f3*f4*f6*f7+k(16)*f1*f3*f5*f6*f7+k(10)*f1*f4*f5*f6*f7+k(10)*f2*f3*f4*f5*f7+k(12)*f2*f3*f4*f6*f7+k(12)*f2*f3*f5*f6*f7+k(10)*f2*f4*f5*f6*f7+k(5)*f1*f2*f3*f7*f8+k(15)//k(2)*f1*f2*f4*f7*f8+k(15)//k(2)*f1*f2*f5*f7*f8+k(33)//k(8)*f1*f2*f6*f7*f8+k(15)//k(2)*f1*f3*f4*f7*f8+k(10)*f1*f3*f5*f7*f8+k(11)//k(2)*f1*f3*f6*f7*f8+k(15)//k(2)*f1*f4*f5*f7*f8+k(11)//k(2)*f1*f4*f6*f7*f8+k(3)*f1*f5*f6*f7*f8+k(5)*f2*f3*f4*f7*f8+k(15)//k(2)*f2*f3*f5*f7*f8+k(33)//k(8)*f2*f3*f6*f7*f8+k(15)//k(2)*f2*f4*f5*f7*f8+k(11)//k(2)*f2*f4*f6*f7*f8+k(3)*f2*f5*f6*f7*f8+k(5)*f3*f4*f5*f7*f8+k(33)//k(8)*f3*f4*f6*f7*f8+k(3)*f3*f5*f6*f7*f8+k(45)//k(8)*f1*f2*f3*f4*f8+k(9)*f1*f2*f3*f5*f8+k(9)*f1*f2*f3*f6*f8+k(81)//k(8)*f1*f2*f4*f5*f8+k(27)//k(2)*f1*f2*f4*f6*f8+k(81)//k(8)*f1*f2*f5*f6*f8+k(9)*f1*f3*f4*f5*f8+k(27)//k(2)*f1*f3*f4*f6*f8+k(27)//k(2)*f1*f3*f5*f6*f8+k(9)*f1*f4*f5*f6*f8+k(45)//k(8)*f2*f3*f4*f5*f8+k(9)*f2*f3*f4*f6*f8+k(81)//k(8)*f2*f3*f5*f6*f8+k(9)*f2*f4*f5*f6*f8+k(45)//k(8)*f3*f4*f5*f6*f8,
        -k(3562995798)+k(5)*f4*f3*f1*f2+k(5)*f2*f3*f4*f5+k(5)*f3*f4*f5*f6+k(8)*f1*f2*f3*f5+k(9)*f1*f2*f4*f5+k(8)*f1*f3*f4*f5+k(5)*f4*f5*f6*f7+k(8)*f1*f2*f3*f6+k(12)*f1*f2*f4*f6+k(9)*f1*f2*f5*f6+k(12)*f1*f3*f4*f6+k(12)*f1*f3*f5*f6+k(8)*f1*f4*f5*f6+k(8)*f2*f3*f4*f6+k(9)*f2*f3*f5*f6+k(8)*f2*f4*f5*f6+k(3)//k(2)*f5*f6*f7*f8+k(8)*f1*f2*f3*f7+k(12)*f1*f2*f4*f7+k(12)*f1*f2*f5*f7+k(9)*f1*f2*f6*f7+k(12)*f1*f3*f4*f7+k(16)*f1*f3*f5*f7+k(12)*f1*f3*f6*f7+k(12)*f1*f4*f5*f7+k(12)*f1*f4*f6*f7+k(8)*f1*f5*f6*f7+k(8)*f2*f3*f4*f7+k(12)*f2*f3*f5*f7+k(9)*f2*f3*f6*f7+k(12)*f2*f4*f5*f7+k(12)*f2*f4*f6*f7+k(8)*f2*f5*f6*f7+k(8)*f3*f4*f5*f7+k(9)*f3*f4*f6*f7+k(8)*f3*f5*f6*f7+k(15)//k(4)*f1*f2*f7*f8+k(5)*f1*f3*f7*f8+k(5)*f1*f4*f7*f8+k(5)*f1*f5*f7*f8+k(11)//k(4)*f1*f6*f7*f8+k(15)//k(4)*f2*f3*f7*f8+k(5)*f2*f4*f7*f8+k(5)*f2*f5*f7*f8+k(11)//k(4)*f2*f6*f7*f8+k(15)//k(4)*f3*f4*f7*f8+k(5)*f3*f5*f7*f8+k(11)//k(4)*f3*f6*f7*f8+k(15)//k(4)*f4*f5*f7*f8+k(11)//k(4)*f4*f6*f7*f8+k(9)//k(2)*f1*f2*f3*f8+k(27)//k(4)*f1*f2*f4*f8+k(27)//k(4)*f1*f2*f5*f8+k(27)//k(4)*f1*f2*f6*f8+k(27)//k(4)*f1*f3*f4*f8+k(9)*f1*f3*f5*f8+k(9)*f1*f3*f6*f8+k(27)//k(4)*f1*f4*f5*f8+k(9)*f1*f4*f6*f8+k(27)//k(4)*f1*f5*f6*f8+k(9)//k(2)*f2*f3*f4*f8+k(27)//k(4)*f2*f3*f5*f8+k(27)//k(4)*f2*f3*f6*f8+k(27)//k(4)*f2*f4*f5*f8+k(9)*f2*f4*f6*f8+k(27)//k(4)*f2*f5*f6*f8+k(9)//k(2)*f3*f4*f5*f8+k(27)//k(4)*f3*f4*f6*f8+k(27)//k(4)*f3*f5*f6*f8+k(9)//k(2)*f4*f5*f6*f8,
        -k(41268600)+k(4)*f2*f3*f4+k(6)*f3*f4*f1+k(4)*f3*f4*f5+k(6)*f1*f2*f4+k(6)*f4*f5*f2+k(6)*f4*f5*f1+k(4)*f4*f5*f6+k(6)*f1*f2*f5+k(8)*f1*f3*f5+k(6)*f2*f3*f5+k(6)*f5*f6*f3+k(6)*f5*f6*f2+k(6)*f5*f6*f1+k(4)*f5*f6*f7+k(6)*f1*f2*f6+k(8)*f1*f3*f6+k(8)*f1*f4*f6+k(6)*f2*f3*f6+k(8)*f2*f4*f6+k(6)*f3*f4*f6+k(6)*f6*f7*f4+k(6)*f6*f7*f3+k(6)*f6*f7*f2+k(6)*f6*f7*f1+k(11)//k(8)*f6*f7*f8+k(6)*f1*f2*f7+k(8)*f1*f3*f7+k(8)*f1*f4*f7+k(8)*f1*f5*f7+k(6)*f2*f3*f7+k(8)*f2*f4*f7+k(8)*f2*f5*f7+k(6)*f3*f4*f7+k(8)*f3*f5*f7+k(6)*f4*f5*f7+k(5)//k(2)*f7*f8*f5+k(5)//k(2)*f7*f8*f4+k(5)//k(2)*f7*f8*f3+k(5)//k(2)*f7*f8*f2+k(5)//k(2)*f7*f8*f1+k(27)//k(8)*f1*f2*f8+k(9)//k(2)*f1*f3*f8+k(9)//k(2)*f1*f4*f8+k(9)//k(2)*f1*f5*f8+k(9)//k(2)*f1*f6*f8+k(27)//k(8)*f2*f3*f8+k(9)//k(2)*f2*f4*f8+k(9)//k(2)*f2*f5*f8+k(9)//k(2)*f2*f6*f8+k(27)//k(8)*f3*f4*f8+k(9)//k(2)*f3*f5*f8+k(9)//k(2)*f3*f6*f8+k(27)//k(8)*f4*f5*f8+k(9)//k(2)*f4*f6*f8+k(27)//k(8)*f5*f6*f8+k(4)*f3*f1*f2,
        -k(257068)+k(3)*f1*f2+k(3)*f2*f3+k(4)*f3*f1+k(3)*f3*f4+k(4)*f4*f2+k(4)*f4*f1+k(3)*f4*f5+k(4)*f5*f3+k(4)*f5*f2+k(4)*f5*f1+k(3)*f5*f6+k(4)*f6*f4+k(4)*f6*f3+k(4)*f6*f2+k(4)*f6*f1+k(3)*f6*f7+k(4)*f7*f5+k(4)*f7*f4+k(4)*f7*f3+k(4)*f7*f2+k(4)*f7*f1+k(5)//k(4)*f7*f8+k(9)//k(4)*f8*f6+k(9)//k(4)*f8*f5+k(9)//k(4)*f8*f4+k(9)//k(4)*f8*f3+k(9)//k(4)*f8*f2+k(9)//k(4)*f8*f1,k(9)//k(8)*f8+k(2)*f7+k(2)*f6+k(2)*f5+k(2)*f4+k(2)*f3+k(2)*f2+k(2)*f1-k(808)
    ]
end

# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion9.ms
function henrion9(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (f1, f2, f3, f4, f5, f6, f7, f8, f9) = np.polynomial_ring(
        k,
        ["f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9"],
        internal_ordering=internal_ordering
    )
    [
        2*f1*f2*f3*f4*f5*f6*f7*f8*f9-22561587455281875,
        9*f8*f7*f6*f5*f4*f3*f1*f2+17//k(9)*f2*f3*f4*f5*f6*f7*f8*f9+77//k(9)*f1*f2*f3*f4*f5*f6*f8*f9+8*f1*f2*f3*f4*f5*f7*f8*f9+65//k(9)*f1*f2*f3*f4*f6*f7*f8*f9+56//k(9)*f1*f2*f3*f5*f6*f7*f8*f9+5*f1*f2*f4*f5*f6*f7*f8*f9+32//k(9)*f1*f3*f4*f5*f6*f7*f8*f9+80//k(9)*f1*f2*f3*f4*f5*f6*f7*f9-10687067741975625,
        -1213257282043500+44//k(3)*f1*f2*f3*f5*f6*f8*f9+32//k(3)*f1*f2*f3*f5*f7*f8*f9+52//k(9)*f1*f2*f3*f6*f7*f8*f9+44//k(3)*f1*f2*f4*f5*f6*f8*f9+12*f1*f2*f4*f5*f7*f8*f9+26//k(3)*f1*f2*f4*f6*f7*f8*f9+14//k(3)*f1*f2*f5*f6*f7*f8*f9+110//k(9)*f1*f3*f4*f5*f6*f8*f9+32//k(3)*f1*f3*f4*f5*f7*f8*f9+26//k(3)*f1*f3*f4*f6*f7*f8*f9+56//k(9)*f1*f3*f5*f6*f7*f8*f9+10//k(3)*f1*f4*f5*f6*f7*f8*f9+22//k(3)*f2*f3*f4*f5*f6*f8*f9+20//k(3)*f2*f3*f4*f5*f7*f8*f9+52//k(9)*f2*f3*f4*f6*f7*f8*f9+14//k(3)*f2*f3*f5*f6*f7*f8*f9+10//k(3)*f2*f4*f5*f6*f7*f8*f9+70//k(9)*f1*f2*f3*f4*f5*f6*f9+40//k(3)*f1*f2*f3*f4*f5*f7*f9+50//k(3)*f1*f2*f3*f4*f6*f7*f9+160//k(9)*f1*f2*f3*f5*f6*f7*f9+50//k(3)*f1*f2*f4*f5*f6*f7*f9+40//k(3)*f1*f3*f4*f5*f6*f7*f9+70//k(9)*f2*f3*f4*f5*f6*f7*f9+8*f7*f6*f5*f4*f3*f1*f2+8*f2*f3*f4*f5*f6*f7*f8+16//k(9)*f3*f4*f5*f6*f7*f8*f9+14*f1*f2*f3*f4*f5*f6*f8+18*f1*f2*f3*f4*f5*f7*f8+20*f1*f2*f3*f4*f6*f7*f8+20*f1*f2*f3*f5*f6*f7*f8+18*f1*f2*f4*f5*f6*f7*f8+14*f1*f3*f4*f5*f6*f7*f8+22//k(3)*f1*f2*f3*f4*f5*f8*f9+110//k(9)*f1*f2*f3*f4*f6*f8*f9+20//k(3)*f1*f2*f3*f4*f7*f8*f9,
        -56374959676500+52//k(9)*f1*f3*f6*f7*f8*f9+100//k(9)*f1*f4*f5*f6*f7*f9+88//k(9)*f1*f4*f5*f6*f8*f9+8*f1*f4*f5*f7*f8*f9+52//k(9)*f1*f4*f6*f7*f8*f9+28//k(9)*f1*f5*f6*f7*f8*f9+20//k(3)*f2*f3*f4*f5*f6*f9+100//k(9)*f2*f3*f4*f5*f7*f9+55//k(9)*f2*f3*f4*f5*f8*f9+40//k(3)*f2*f3*f4*f6*f7*f9+88//k(9)*f2*f3*f4*f6*f8*f9+16//k(3)*f2*f3*f4*f7*f8*f9+40//k(3)*f2*f3*f5*f6*f7*f9+11*f2*f3*f5*f6*f8*f9+8*f2*f3*f5*f7*f8*f9+13//k(3)*f2*f3*f6*f7*f8*f9+100//k(9)*f2*f4*f5*f6*f7*f9+88//k(9)*f2*f4*f5*f6*f8*f9+8*f2*f4*f5*f7*f8*f9+52//k(9)*f2*f4*f6*f7*f8*f9+28//k(9)*f2*f5*f6*f7*f8*f9+20//k(3)*f3*f4*f5*f6*f7*f9+55//k(9)*f3*f4*f5*f6*f8*f9+16//k(3)*f3*f4*f5*f7*f8*f9+13//k(3)*f3*f4*f6*f7*f8*f9+28//k(9)*f3*f5*f6*f7*f8*f9+20*f1*f3*f4*f5*f6*f8+24*f1*f3*f4*f5*f7*f8+24*f1*f3*f4*f6*f7*f8+20*f1*f3*f5*f6*f7*f8+12*f1*f4*f5*f6*f7*f8+12*f2*f3*f4*f5*f6*f8+15*f2*f3*f4*f5*f7*f8+16*f2*f3*f4*f6*f7*f8+15*f2*f3*f5*f6*f7*f8+12*f2*f4*f5*f6*f7*f8+20//k(3)*f1*f2*f3*f4*f5*f9+100//k(9)*f1*f2*f3*f4*f6*f9+100//k(9)*f1*f2*f3*f4*f7*f9+55//k(9)*f1*f2*f3*f4*f8*f9+40//k(3)*f1*f2*f3*f5*f6*f9+160//k(9)*f1*f2*f3*f5*f7*f9+88//k(9)*f1*f2*f3*f5*f8*f9+40//k(3)*f1*f2*f3*f6*f7*f9+88//k(9)*f1*f2*f3*f6*f8*f9+16//k(3)*f1*f2*f3*f7*f8*f9+40//k(3)*f1*f2*f4*f5*f6*f9+20*f1*f2*f4*f5*f7*f9+11*f1*f2*f4*f5*f8*f9+20*f1*f2*f4*f6*f7*f9+44//k(3)*f1*f2*f4*f6*f8*f9+8*f1*f2*f4*f7*f8*f9+40//k(3)*f1*f2*f5*f6*f7*f9+11*f1*f2*f5*f6*f8*f9+8*f1*f2*f5*f7*f8*f9+13//k(3)*f1*f2*f6*f7*f8*f9+100//k(9)*f1*f3*f4*f5*f6*f9+160//k(9)*f1*f3*f4*f5*f7*f9+88//k(9)*f1*f3*f4*f5*f8*f9+20*f1*f3*f4*f6*f7*f9+44//k(3)*f1*f3*f4*f6*f8*f9+8*f1*f3*f4*f7*f8*f9+160//k(9)*f1*f3*f5*f6*f7*f9+44//k(3)*f1*f3*f5*f6*f8*f9+32//k(3)*f1*f3*f5*f7*f8*f9+7*f3*f4*f5*f6*f7*f8+12*f1*f2*f3*f4*f5*f7+15*f1*f2*f3*f4*f6*f7+16*f1*f2*f3*f5*f6*f7+15*f1*f2*f4*f5*f6*f7+12*f1*f3*f4*f5*f6*f7+5//k(3)*f4*f5*f6*f7*f8*f9+12*f1*f2*f3*f4*f5*f8+20*f1*f2*f3*f4*f6*f8+15*f1*f2*f3*f4*f7*f8+24*f1*f2*f3*f5*f6*f8+24*f1*f2*f3*f5*f7*f8+16*f1*f2*f3*f6*f7*f8+24*f1*f2*f4*f5*f6*f8+27*f1*f2*f4*f5*f7*f8+24*f1*f2*f4*f6*f7*f8+15*f1*f2*f5*f6*f7*f8+7*f6*f5*f4*f3*f1*f2+7*f2*f3*f4*f5*f6*f7,
        -1314069041754+6*f5*f4*f3*f1*f2+6*f2*f3*f4*f5*f6+6*f3*f4*f5*f6*f7+10*f1*f2*f3*f4*f6+12*f1*f2*f3*f5*f6+12*f1*f2*f4*f5*f6+10*f1*f3*f4*f5*f6+6*f4*f5*f6*f7*f8+10*f1*f2*f3*f4*f7+16*f1*f2*f3*f5*f7+12*f1*f2*f3*f6*f7+18*f1*f2*f4*f5*f7+18*f1*f2*f4*f6*f7+12*f1*f2*f5*f6*f7+16*f1*f3*f4*f5*f7+18*f1*f3*f4*f6*f7+16*f1*f3*f5*f6*f7+10*f1*f4*f5*f6*f7+10*f2*f3*f4*f5*f7+12*f2*f3*f4*f6*f7+12*f2*f3*f5*f6*f7+10*f2*f4*f5*f6*f7+14//k(9)*f5*f6*f7*f8*f9+10*f1*f2*f3*f4*f8+16*f1*f2*f3*f5*f8+16*f1*f2*f3*f6*f8+12*f1*f2*f3*f7*f8+18*f1*f2*f4*f5*f8+24*f1*f2*f4*f6*f8+18*f1*f2*f4*f7*f8+18*f1*f2*f5*f6*f8+18*f1*f2*f5*f7*f8+12*f1*f2*f6*f7*f8+16*f1*f3*f4*f5*f8+24*f1*f3*f4*f6*f8+18*f1*f3*f4*f7*f8+24*f1*f3*f5*f6*f8+24*f1*f3*f5*f7*f8+16*f1*f3*f6*f7*f8+16*f1*f4*f5*f6*f8+18*f1*f4*f5*f7*f8+16*f1*f4*f6*f7*f8+10*f1*f5*f6*f7*f8+10*f2*f3*f4*f5*f8+16*f2*f3*f4*f6*f8+12*f2*f3*f4*f7*f8+18*f2*f3*f5*f6*f8+18*f2*f3*f5*f7*f8+12*f2*f3*f6*f7*f8+16*f2*f4*f5*f6*f8+18*f2*f4*f5*f7*f8+16*f2*f4*f6*f7*f8+10*f2*f5*f6*f7*f8+10*f3*f4*f5*f6*f8+12*f3*f4*f5*f7*f8+12*f3*f4*f6*f7*f8+10*f3*f5*f6*f7*f8+50//k(9)*f1*f2*f3*f4*f9+80//k(9)*f1*f2*f3*f5*f9+80//k(9)*f1*f2*f3*f6*f9+80//k(9)*f1*f2*f3*f7*f9+44//k(9)*f1*f2*f3*f8*f9+10*f1*f2*f4*f5*f9+40//k(3)*f1*f2*f4*f6*f9+40//k(3)*f1*f2*f4*f7*f9+22//k(3)*f1*f2*f4*f8*f9+10*f1*f2*f5*f6*f9+40//k(3)*f1*f2*f5*f7*f9+22//k(3)*f1*f2*f5*f8*f9+10*f1*f2*f6*f7*f9+22//k(3)*f1*f2*f6*f8*f9+4*f1*f2*f7*f8*f9+80//k(9)*f1*f3*f4*f5*f9+40//k(3)*f1*f3*f4*f6*f9+40//k(3)*f1*f3*f4*f7*f9+22//k(3)*f1*f3*f4*f8*f9+40//k(3)*f1*f3*f5*f6*f9+160//k(9)*f1*f3*f5*f7*f9+88//k(9)*f1*f3*f5*f8*f9+40//k(3)*f1*f3*f6*f7*f9+88//k(9)*f1*f3*f6*f8*f9+16//k(3)*f1*f3*f7*f8*f9+80//k(9)*f1*f4*f5*f6*f9+40//k(3)*f1*f4*f5*f7*f9+22//k(3)*f1*f4*f5*f8*f9+40//k(3)*f1*f4*f6*f7*f9+88//k(9)*f1*f4*f6*f8*f9+16//k(3)*f1*f4*f7*f8*f9+80//k(9)*f1*f5*f6*f7*f9+22//k(3)*f1*f5*f6*f8*f9+16//k(3)*f1*f5*f7*f8*f9+26//k(9)*f1*f6*f7*f8*f9+50//k(9)*f2*f3*f4*f5*f9+80//k(9)*f2*f3*f4*f6*f9+80//k(9)*f2*f3*f4*f7*f9+44//k(9)*f2*f3*f4*f8*f9+10*f2*f3*f5*f6*f9+40//k(3)*f2*f3*f5*f7*f9+22//k(3)*f2*f3*f5*f8*f9+10*f2*f3*f6*f7*f9+22//k(3)*f2*f3*f6*f8*f9+4*f2*f3*f7*f8*f9+80//k(9)*f2*f4*f5*f6*f9+40//k(3)*f2*f4*f5*f7*f9+22//k(3)*f2*f4*f5*f8*f9+40//k(3)*f2*f4*f6*f7*f9+88//k(9)*f2*f4*f6*f8*f9+16//k(3)*f2*f4*f7*f8*f9+80//k(9)*f2*f5*f6*f7*f9+22//k(3)*f2*f5*f6*f8*f9+16//k(3)*f2*f5*f7*f8*f9+26//k(9)*f2*f6*f7*f8*f9+50//k(9)*f3*f4*f5*f6*f9+80//k(9)*f3*f4*f5*f7*f9+44//k(9)*f3*f4*f5*f8*f9+10*f3*f4*f6*f7*f9+22//k(3)*f3*f4*f6*f8*f9+4*f3*f4*f7*f8*f9+80//k(9)*f3*f5*f6*f7*f9+22//k(3)*f3*f5*f6*f8*f9+16//k(3)*f3*f5*f7*f8*f9+26//k(9)*f3*f6*f7*f8*f9+50//k(9)*f4*f5*f6*f7*f9+44//k(9)*f4*f5*f6*f8*f9+4*f4*f5*f7*f8*f9+26//k(9)*f4*f6*f7*f8*f9,
        -16892753598+5*f4*f3*f1*f2+5*f2*f3*f4*f5+5*f3*f4*f5*f6+8*f1*f2*f3*f5+9*f1*f2*f4*f5+8*f1*f3*f4*f5+5*f4*f5*f6*f7+8*f1*f2*f3*f6+12*f1*f2*f4*f6+9*f1*f2*f5*f6+12*f1*f3*f4*f6+12*f1*f3*f5*f6+8*f1*f4*f5*f6+8*f2*f3*f4*f6+9*f2*f3*f5*f6+8*f2*f4*f5*f6+5*f5*f6*f7*f8+8*f1*f2*f3*f7+12*f1*f2*f4*f7+12*f1*f2*f5*f7+9*f1*f2*f6*f7+12*f1*f3*f4*f7+16*f1*f3*f5*f7+12*f1*f3*f6*f7+12*f1*f4*f5*f7+12*f1*f4*f6*f7+8*f1*f5*f6*f7+8*f2*f3*f4*f7+12*f2*f3*f5*f7+9*f2*f3*f6*f7+12*f2*f4*f5*f7+12*f2*f4*f6*f7+8*f2*f5*f6*f7+8*f3*f4*f5*f7+9*f3*f4*f6*f7+8*f3*f5*f6*f7+13//k(9)*f6*f7*f8*f9+8*f1*f2*f3*f8+12*f1*f2*f4*f8+12*f1*f2*f5*f8+12*f1*f2*f6*f8+9*f1*f2*f7*f8+12*f1*f3*f4*f8+16*f1*f3*f5*f8+16*f1*f3*f6*f8+12*f1*f3*f7*f8+12*f1*f4*f5*f8+16*f1*f4*f6*f8+12*f1*f4*f7*f8+12*f1*f5*f6*f8+12*f1*f5*f7*f8+8*f1*f6*f7*f8+8*f2*f3*f4*f8+12*f2*f3*f5*f8+12*f2*f3*f6*f8+9*f2*f3*f7*f8+12*f2*f4*f5*f8+16*f2*f4*f6*f8+12*f2*f4*f7*f8+12*f2*f5*f6*f8+12*f2*f5*f7*f8+8*f2*f6*f7*f8+8*f3*f4*f5*f8+12*f3*f4*f6*f8+9*f3*f4*f7*f8+12*f3*f5*f6*f8+12*f3*f5*f7*f8+8*f3*f6*f7*f8+8*f4*f5*f6*f8+9*f4*f5*f7*f8+8*f4*f6*f7*f8+40//k(9)*f1*f2*f3*f9+20//k(3)*f1*f2*f4*f9+20//k(3)*f1*f2*f5*f9+20//k(3)*f1*f2*f6*f9+20//k(3)*f1*f2*f7*f9+11//k(3)*f1*f2*f8*f9+20//k(3)*f1*f3*f4*f9+80//k(9)*f1*f3*f5*f9+80//k(9)*f1*f3*f6*f9+80//k(9)*f1*f3*f7*f9+44//k(9)*f1*f3*f8*f9+20//k(3)*f1*f4*f5*f9+80//k(9)*f1*f4*f6*f9+80//k(9)*f1*f4*f7*f9+44//k(9)*f1*f4*f8*f9+20//k(3)*f1*f5*f6*f9+80//k(9)*f1*f5*f7*f9+44//k(9)*f1*f5*f8*f9+20//k(3)*f1*f6*f7*f9+44//k(9)*f1*f6*f8*f9+8//k(3)*f1*f7*f8*f9+40//k(9)*f2*f3*f4*f9+20//k(3)*f2*f3*f5*f9+20//k(3)*f2*f3*f6*f9+20//k(3)*f2*f3*f7*f9+11//k(3)*f2*f3*f8*f9+20//k(3)*f2*f4*f5*f9+80//k(9)*f2*f4*f6*f9+80//k(9)*f2*f4*f7*f9+44//k(9)*f2*f4*f8*f9+20//k(3)*f2*f5*f6*f9+80//k(9)*f2*f5*f7*f9+44//k(9)*f2*f5*f8*f9+20//k(3)*f2*f6*f7*f9+44//k(9)*f2*f6*f8*f9+8//k(3)*f2*f7*f8*f9+40//k(9)*f3*f4*f5*f9+20//k(3)*f3*f4*f6*f9+20//k(3)*f3*f4*f7*f9+11//k(3)*f3*f4*f8*f9+20//k(3)*f3*f5*f6*f9+80//k(9)*f3*f5*f7*f9+44//k(9)*f3*f5*f8*f9+20//k(3)*f3*f6*f7*f9+44//k(9)*f3*f6*f8*f9+8//k(3)*f3*f7*f8*f9+40//k(9)*f4*f5*f6*f9+20//k(3)*f4*f5*f7*f9+11//k(3)*f4*f5*f8*f9+20//k(3)*f4*f6*f7*f9+44//k(9)*f4*f6*f8*f9+8//k(3)*f4*f7*f8*f9+40//k(9)*f5*f6*f7*f9+11//k(3)*f5*f6*f8*f9+8//k(3)*f5*f7*f8*f9,
        -124301564+4*f3*f1*f2+4*f2*f3*f4+6*f3*f4*f1+4*f3*f4*f5+6*f1*f2*f4+6*f4*f5*f2+6*f4*f5*f1+4*f4*f5*f6+6*f1*f2*f5+8*f1*f3*f5+6*f2*f3*f5+6*f5*f6*f3+6*f5*f6*f2+6*f5*f6*f1+4*f5*f6*f7+6*f1*f2*f6+8*f1*f3*f6+8*f1*f4*f6+6*f2*f3*f6+8*f2*f4*f6+6*f3*f4*f6+6*f6*f7*f4+6*f6*f7*f3+6*f6*f7*f2+6*f6*f7*f1+4*f6*f7*f8+6*f1*f2*f7+8*f1*f3*f7+8*f1*f4*f7+8*f1*f5*f7+6*f2*f3*f7+8*f2*f4*f7+8*f2*f5*f7+6*f3*f4*f7+8*f3*f5*f7+6*f4*f5*f7+6*f7*f8*f5+6*f7*f8*f4+6*f7*f8*f3+6*f7*f8*f2+6*f7*f8*f1+4//k(3)*f7*f8*f9+6*f1*f2*f8+8*f1*f3*f8+8*f1*f4*f8+8*f1*f5*f8+8*f1*f6*f8+6*f2*f3*f8+8*f2*f4*f8+8*f2*f5*f8+8*f2*f6*f8+6*f3*f4*f8+8*f3*f5*f8+8*f3*f6*f8+6*f4*f5*f8+8*f4*f6*f8+6*f5*f6*f8+22//k(9)*f8*f9*f6+22//k(9)*f8*f9*f5+22//k(9)*f8*f9*f4+22//k(9)*f8*f9*f3+22//k(9)*f8*f9*f2+22//k(9)*f8*f9*f1+10//k(3)*f1*f2*f9+40//k(9)*f1*f3*f9+40//k(9)*f1*f4*f9+40//k(9)*f1*f5*f9+40//k(9)*f1*f6*f9+40//k(9)*f1*f7*f9+10//k(3)*f2*f3*f9+40//k(9)*f2*f4*f9+40//k(9)*f2*f5*f9+40//k(9)*f2*f6*f9+40//k(9)*f2*f7*f9+10//k(3)*f3*f4*f9+40//k(9)*f3*f5*f9+40//k(9)*f3*f6*f9+40//k(9)*f3*f7*f9+10//k(3)*f4*f5*f9+40//k(9)*f4*f6*f9+40//k(9)*f4*f7*f9+10//k(3)*f5*f6*f9+40//k(9)*f5*f7*f9+10//k(3)*f6*f7*f9,
        -518052+3*f1*f2+3*f2*f3+4*f3*f1+3*f3*f4+4*f4*f2+4*f4*f1+3*f4*f5+4*f5*f3+4*f5*f2+4*f5*f1+3*f5*f6+4*f6*f4+4*f6*f3+4*f6*f2+4*f6*f1+3*f6*f7+4*f7*f5+4*f7*f4+4*f7*f3+4*f7*f2+4*f7*f1+3*f7*f8+4*f8*f6+4*f8*f5+4*f8*f4+4*f8*f3+4*f8*f2+4*f8*f1+11//k(9)*f8*f9+20//k(9)*f9*f7+20//k(9)*f9*f6+20//k(9)*f9*f5+20//k(9)*f9*f4+20//k(9)*f9*f3+20//k(9)*f9*f2+20//k(9)*f9*f1,
        10//k(9)*f9+2*f8+2*f7+2*f6+2*f5+2*f4+2*f3+2*f2+2*f1-1131
    ]
end

###
# Eco

function econ(n; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
	R, x = polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=internal_ordering)
	vcat(
		[
			x[n]*(x[k] + sum(x[1:n-k-1] .* x[1+k:n-1]; init=0)) - k
			for k in 1:n-1
		],
		sum(x[1:n-1]; init=0) + 1
    )
end

function eco5(; np=AbstractAlgebra, k=np.GF(2^31 - 1), internal_ordering=:degrevlex)
    _, (x1, x2, x3, x4, x5) =
        np.polynomial_ring(k, ["x$i" for i in 1:5], internal_ordering=internal_ordering)
    [
        (x1 + x1 * x2 + x2 * x3 + x3 * x4) * x5 - 1,
        (x2 + x1 * x3 + x2 * x4) * x5 - 2,
        (x3 + x1 * x4) * x5 - 3,
        x4 * x5 - 4,
        x1 + x2 + x3 + x4 + 1
    ]
end

function eco7(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x1, x2, x3, x4, x5, x6, x7) =
        np.polynomial_ring(k, ["x$i" for i in 1:7], internal_ordering=internal_ordering)
    [
        (x1 + x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x5 * x6) * x7 - 1,
        (x2 + x1 * x3 + x2 * x4 + x3 * x5 + x4 * x6) * x7 - 2,
        (x3 + x1 * x4 + x2 * x5 + x3 * x6) * x7 - 3,
        (x4 + x1 * x5 + x2 * x6) * x7 - 4,
        (x5 + x1 * x6) * x7 - 5,
        x6 * x7 - 6,
        x1 + x2 + x3 + x4 + x5 + x6 + 1
    ]
end

function eco10(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9) =
        np.polynomial_ring(k, ["x$i" for i in 0:9], internal_ordering=internal_ordering)
    [
        x0 * x1 * x9 + x1 * x2 * x9 + x2 * x3 * x9 + x3 * x4 * x9 + x4 * x5 * x9 + x5 * x6 * x9 + x6 * x7 * x9 + x7 * x8 * x9 + x0 * x9 - 1, 
        x0 * x2 * x9 + x1 * x3 * x9 + x2 * x4 * x9 + x3 * x5 * x9 + x4 * x6 * x9 + x5 * x7 * x9 + x6 * x8 * x9 + x1 * x9 - 2, 
        x0 * x3 * x9 + x1 * x4 * x9 + x2 * x5 * x9 + x3 * x6 * x9 + x4 * x7 * x9 + x5 * x8 * x9 + x2 * x9 - 3,
        x0 * x4 * x9 + x1 * x5 * x9 + x2 * x6 * x9 + x3 * x7 * x9 + x4 * x8 * x9 + x3 * x9 - 4,
        x0 * x5 * x9 + x1 * x6 * x9 + x2 * x7 * x9 + x3 * x8 * x9 + x4 * x9 - 5,
        x0 * x6 * x9 + x1 * x7 * x9 + x2 * x8 * x9 + x5 * x9 - 6,
        x0 * x7 * x9 + x1 * x8 * x9 + x6 * x9 - 7,
        x0 * x8 * x9 + x7 * x9 - 8,
        x8 * x9 - 9,
        x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + 1
    ]
end

# Critical pair degree sequence: 3,3,3,3,4,3,4,5,5,5,6,7,7,7,7,7,7,7,8,9
function eco11(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) =
        np.polynomial_ring(k, ["x$i" for i in 0:10], internal_ordering=internal_ordering)
    [
        x0 * x1 * x10 + x1 * x2 * x10 + x2 * x3 * x10 + x3 * x4 * x10 + x4 * x5 * x10 + x5 * x6 * x10 + x6 * x7 * x10 + x7 * x8 * x10 + x8 * x9 * x10 + x0 * x10 - 1, 
        x0 * x2 * x10 + x1 * x3 * x10 + x2 * x4 * x10 + x3 * x5 * x10 + x4 * x6 * x10 + x5 * x7 * x10 + x6 * x8 * x10 + x7 * x9 * x10 + x1 * x10 - 2,
        x0 * x3 * x10 + x1 * x4 * x10 + x2 * x5 * x10 + x3 * x6 * x10 + x4 * x7 * x10 + x5 * x8 * x10 + x6 * x9 * x10 + x2 * x10 - 3,
        x0 * x4 * x10 + x1 * x5 * x10 + x2 * x6 * x10 + x3 * x7 * x10 + x4 * x8 * x10 + x5 * x9 * x10 + x3 * x10 - 4,
        x0 * x5 * x10 + x1 * x6 * x10 + x2 * x7 * x10 + x3 * x8 * x10 + x4 * x9 * x10 + x4 * x10 - 5,
        x0 * x6 * x10 + x1 * x7 * x10 + x2 * x8 * x10 + x3 * x9 * x10 + x5 * x10 - 6,
        x0 * x7 * x10 + x1 * x8 * x10 + x2 * x9 * x10 + x6 * x10 - 7,
        x0 * x8 * x10 + x1 * x9 * x10 + x7 * x10 - 8,
        x0 * x9 * x10 + x8 * x10 - 9,
        x9 * x10 - 10,
        x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 1
    ]
end

function eco12(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) =
        np.polynomial_ring(k, ["x$i" for i in 0:11], internal_ordering=internal_ordering)
    [ 
        x0 * x1 * x11 + x1 * x2 * x11 + x2 * x3 * x11 + x3 * x4 * x11 + x4 * x5 * x11 + x5 * x6 * x11 + x6 * x7 * x11 + x7 * x8 * x11 + x8 * x9 * x11 + x9 * x10 * x11 + x0 * x11 - 1,
        x0 * x2 * x11 + x1 * x3 * x11 + x2 * x4 * x11 + x3 * x5 * x11 + x4 * x6 * x11 + x5 * x7 * x11 + x6 * x8 * x11 + x7 * x9 * x11 + x8 * x10 * x11 + x1 * x11 - 2,
        x0 * x3 * x11 + x1 * x4 * x11 + x2 * x5 * x11 + x3 * x6 * x11 + x4 * x7 * x11 + x5 * x8 * x11 + x6 * x9 * x11 + x7 * x10 * x11 + x2 * x11 - 3,
        x0 * x4 * x11 + x1 * x5 * x11 + x2 * x6 * x11 + x3 * x7 * x11 + x4 * x8 * x11 + x5 * x9 * x11 + x6 * x10 * x11 + x3 * x11 - 4,
        x0 * x5 * x11 + x1 * x6 * x11 + x2 * x7 * x11 + x3 * x8 * x11 + x4 * x9 * x11 + x5 * x10 * x11 + x4 * x11 - 5, 
        x0 * x6 * x11 + x1 * x7 * x11 + x2 * x8 * x11 + x3 * x9 * x11 + x4 * x10 * x11 + x5 * x11 - 6,
        x0 * x7 * x11 + x1 * x8 * x11 + x2 * x9 * x11 + x3 * x10 * x11 + x6 * x11 - 7,
        x0 * x8 * x11 + x1 * x9 * x11 + x2 * x10 * x11 + x7 * x11 - 8,
        x0 * x9 * x11 + x1 * x10 * x11 + x8 * x11 - 9,
        x0 * x10 * x11 + x9 * x11 - 10,
        x10 * x11 - 11,
        x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + 1
    ]
end

function eco13(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12) =
        np.polynomial_ring(k, ["x$i" for i in 0:12], internal_ordering=internal_ordering)
    [
        x0 * x1 * x12 + x1 * x2 * x12 + x2 * x3 * x12 + x3 * x4 * x12 + x4 * x5 * x12 + x5 * x6 * x12 + x6 * x7 * x12 + x7 * x8 * x12 + x8 * x9 * x12 + x9 * x10 * x12 + x10 * x11 * x12 + x0 * x12 - 1,
        x0 * x2 * x12 + x1 * x3 * x12 + x2 * x4 * x12 + x3 * x5 * x12 + x4 * x6 * x12 + x5 * x7 * x12 + x6 * x8 * x12 + x7 * x9 * x12 + x8 * x10 * x12 + x9 * x11 * x12 + x1 * x12 - 2,
        x0 * x3 * x12 + x1 * x4 * x12 + x2 * x5 * x12 + x3 * x6 * x12 + x4 * x7 * x12 + x5 * x8 * x12 + x6 * x9 * x12 + x7 * x10 * x12 + x8 * x11 * x12 + x2 * x12 - 3,
        x0 * x4 * x12 + x1 * x5 * x12 + x2 * x6 * x12 + x3 * x7 * x12 + x4 * x8 * x12 + x5 * x9 * x12 + x6 * x10 * x12 + x7 * x11 * x12 + x3 * x12 - 4,
        x0 * x5 * x12 + x1 * x6 * x12 + x2 * x7 * x12 + x3 * x8 * x12 + x4 * x9 * x12 + x5 * x10 * x12 + x6 * x11 * x12 + x4 * x12 - 5,
        x0 * x6 * x12 + x1 * x7 * x12 + x2 * x8 * x12 + x3 * x9 * x12 + x4 * x10 * x12 + x5 * x11 * x12 + x5 * x12 - 6,
        x0 * x7 * x12 + x1 * x8 * x12 + x2 * x9 * x12 + x3 * x10 * x12 + x4 * x11 * x12 + x6 * x12 - 7,
        x0 * x8 * x12 + x1 * x9 * x12 + x2 * x10 * x12 + x3 * x11 * x12 + x7 * x12 - 8,
        x0 * x9 * x12 + x1 * x10 * x12 + x2 * x11 * x12 + x8 * x12 - 9,
        x0 * x10 * x12 + x1 * x11 * x12 + x9 * x12 - 10,
        x0 * x11 * x12 + x10 * x12 - 11,
        x11 * x12 - 12,
        x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + 1
    ]
end

function eco14(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12, x13) = np.polynomial_ring(k, ["x$i" for i in 0:13], internal_ordering=internal_ordering)
    [
        x0*x1*x13+x1*x2*x13+x2*x3*x13+x3*x4*x13+x4*x5*x13+x5*x6*x13+x6*x7*x13+x7*x8*x13+x8*x9*x13+x9*x10*x13+x10*x11*x13+x11*x12*x13+x0*x13-1,
        x0*x2*x13+x1*x3*x13+x2*x4*x13+x3*x5*x13+x4*x6*x13+x5*x7*x13+x6*x8*x13+x7*x9*x13+x8*x10*x13+x9*x11*x13+x10*x12*x13+x1*x13-2,
        x0*x3*x13+x1*x4*x13+x2*x5*x13+x3*x6*x13+x4*x7*x13+x5*x8*x13+x6*x9*x13+x7*x10*x13+x8*x11*x13+x9*x12*x13+x2*x13-3,
        x0*x4*x13+x1*x5*x13+x2*x6*x13+x3*x7*x13+x4*x8*x13+x5*x9*x13+x6*x10*x13+x7*x11*x13+x8*x12*x13+x3*x13-4,
        x0*x5*x13+x1*x6*x13+x2*x7*x13+x3*x8*x13+x4*x9*x13+x5*x10*x13+x6*x11*x13+x7*x12*x13+x4*x13-5,
        x0*x6*x13+x1*x7*x13+x2*x8*x13+x3*x9*x13+x4*x10*x13+x5*x11*x13+x6*x12*x13+x5*x13-6,
        x0*x7*x13+x1*x8*x13+x2*x9*x13+x3*x10*x13+x4*x11*x13+x5*x12*x13+x6*x13-7,
        x0*x8*x13+x1*x9*x13+x2*x10*x13+x3*x11*x13+x4*x12*x13+x7*x13-8,
        x0*x9*x13+x1*x10*x13+x2*x11*x13+x3*x12*x13+x8*x13-9,
        x0*x10*x13+x1*x11*x13+x2*x12*x13+x9*x13-10,
        x0*x11*x13+x1*x12*x13+x10*x13-11,
        x0*x12*x13+x11*x13-12,
        x12*x13-13,
        x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+1
    ]
end

function eco15(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14) = np.polynomial_ring(k, ["x$i" for i in 0:14], internal_ordering=internal_ordering)
    [
        x0*x1*x14+x1*x2*x14+x2*x3*x14+x3*x4*x14+x4*x5*x14+x5*x6*x14+x6*x7*x14+x7*x8*x14+x8*x9*x14+x9*x10*x14+x10*x11*x14+x11*x12*x14+x12*x13*x14+x0*x14-1,
        x0*x2*x14+x1*x3*x14+x2*x4*x14+x3*x5*x14+x4*x6*x14+x5*x7*x14+x6*x8*x14+x7*x9*x14+x8*x10*x14+x9*x11*x14+x10*x12*x14+x11*x13*x14+x1*x14-2,
        x0*x3*x14+x1*x4*x14+x2*x5*x14+x3*x6*x14+x4*x7*x14+x5*x8*x14+x6*x9*x14+x7*x10*x14+x8*x11*x14+x9*x12*x14+x10*x13*x14+x2*x14-3,
        x0*x4*x14+x1*x5*x14+x2*x6*x14+x3*x7*x14+x4*x8*x14+x5*x9*x14+x6*x10*x14+x7*x11*x14+x8*x12*x14+x9*x13*x14+x3*x14-4,
        x0*x5*x14+x1*x6*x14+x2*x7*x14+x3*x8*x14+x4*x9*x14+x5*x10*x14+x6*x11*x14+x7*x12*x14+x8*x13*x14+x4*x14-5,
        x0*x6*x14+x1*x7*x14+x2*x8*x14+x3*x9*x14+x4*x10*x14+x5*x11*x14+x6*x12*x14+x7*x13*x14+x5*x14-6,
        x0*x7*x14+x1*x8*x14+x2*x9*x14+x3*x10*x14+x4*x11*x14+x5*x12*x14+x6*x13*x14+x6*x14-7,
        x0*x8*x14+x1*x9*x14+x2*x10*x14+x3*x11*x14+x4*x12*x14+x5*x13*x14+x7*x14-8,
        x0*x9*x14+x1*x10*x14+x2*x11*x14+x3*x12*x14+x4*x13*x14+x8*x14-9,
        x0*x10*x14+x1*x11*x14+x2*x12*x14+x3*x13*x14+x9*x14-10,
        x0*x11*x14+x1*x12*x14+x2*x13*x14+x10*x14-11,
        x0*x12*x14+x1*x13*x14+x11*x14-12,
        x0*x13*x14+x12*x14-13,
        x13*x14-14,
        x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+1
    ]
end

# Source:
# https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/master/src/rps10.jl
function rps10(; tol=0, np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (g1, g2, g3, p1, p2, p3, q0, q1, q2, q3) = np.polynomial_ring(k, ["g1", "g2", "g3", "p1", "p2", "p3", "q0", "q1", "q2", "q3",], internal_ordering=internal_ordering)
    equations = [
    -rationalize(BigInt, 0.1279703687075118, tol=tol)*g1^2 - rationalize(BigInt, 0.48596123125526264, tol=tol)*g1*g2 + rationalize(BigInt, 0.30699556370717496, tol=tol)*g2^2 + rationalize(BigInt, 0.3778977698527674, tol=tol)*g1*g3 - rationalize(BigInt, 0.23404544076569642, tol=tol)*g2*g3 + rationalize(BigInt, 0.01563626178508072, tol=tol)*g3^2 + rationalize(BigInt, 0.327228678790004, tol=tol)*g1^2*p1 + rationalize(BigInt, 0.8426829275672494, tol=tol)*g1*g2*p1 + rationalize(BigInt, 0.6075645757034159, tol=tol)*g2^2*p1 - rationalize(BigInt, 1.1371405598667543, tol=tol)*g1*g3*p1 + rationalize(BigInt, 0.229293271620915, tol=tol)*g2*g3*p1 - rationalize(BigInt, 0.21948911177437957, tol=tol)*g3^2*p1 - rationalize(BigInt, 0.2075154964282774, tol=tol)*g1^2*p1^2 - rationalize(BigInt, 0.37702968479068544, tol=tol)*g1*g2*p1^2 - rationalize(BigInt, 0.16688906819159421, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 0.7986954318323025, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 0.866826144775651, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.37440456461987165, tol=tol)*g3^2*p1^2 + rationalize(BigInt, 1.5614616440131446, tol=tol)*g1^2*p2 - rationalize(BigInt, 1.7388380675822595, tol=tol)*g1*g2*p2 + rationalize(BigInt, 0.06790915713070725, tol=tol)*g2^2*p2 - rationalize(BigInt, 0.4309121044684771, tol=tol)*g1*g3*p2 + rationalize(BigInt, 0.9086272006283425, tol=tol)*g2*g3*p2 - rationalize(BigInt, 0.2764931751394603, tol=tol)*g3^2*p2 - rationalize(BigInt, 1.8163349832174116, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 0.9167144057621401, tol=tol)*g1*g2*p1*p2 + rationalize(BigInt, 1.0203368504488892, tol=tol)*g2^2*p1*p2 - rationalize(BigInt, 0.23194646823111892, tol=tol)*g1*g3*p1*p2 + rationalize(BigInt, 0.539670777307627, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 0.7959981327685224, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.08717268867521591, tol=tol)*g1^2*p2^2 + rationalize(BigInt, 0.9504154644263471, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.48206756571420756, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 1.065062423127697, tol=tol)*g1*g3*p2^2 + rationalize(BigInt, 0.1209952909274163, tol=tol)*g2*g3*p2^2 + rationalize(BigInt, 0.3948948770389917, tol=tol)*g3^2*p2^2 + rationalize(BigInt, 0.289766299873838, tol=tol)*g1^2*p3 - rationalize(BigInt, 1.2778927965251532, tol=tol)*g1*g2*p3 + rationalize(BigInt, 0.9087896778886251, tol=tol)*g2^2*p3 - rationalize(BigInt, 0.5812612591154215, tol=tol)*g1*g3*p3 - rationalize(BigInt, 0.7595904624983555, tol=tol)*g2*g3*p3 + rationalize(BigInt, 0.5084892760496751, tol=tol)*g3^2*p3 - rationalize(BigInt, 0.3268802641947883, tol=tol)*g1^2*p1*p3 + rationalize(BigInt, 0.657630238424344, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 1.1093919363972093, tol=tol)*g2^2*p1*p3 + rationalize(BigInt, 0.4551393419480071, tol=tol)*g1*g3*p1*p3 + rationalize(BigInt, 1.8553852513069364, tol=tol)*g2*g3*p1*p3 - rationalize(BigInt, 0.7825116722024211, tol=tol)*g3^2*p1*p3 + rationalize(BigInt, 0.5810469298461638, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 1.1557382363783264, tol=tol)*g1*g2*p2*p3 - rationalize(BigInt, 0.11367961187637783, tol=tol)*g2^2*p2*p3 + rationalize(BigInt, 1.7077140933509898, tol=tol)*g1*g3*p2*p3 - rationalize(BigInt, 0.36547942767108677, tol=tol)*g2*g3*p2*p3 - rationalize(BigInt, 0.4673673179697859, tol=tol)*g3^2*p2*p3 + rationalize(BigInt, 0.12034280775306151, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 0.5733857796356615, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.6489566339058018, tol=tol)*g2^2*p3^2 + rationalize(BigInt, 0.2663669912953945, tol=tol)*g1*g3*p3^2 - rationalize(BigInt, 0.9878214357030672, tol=tol)*g2*g3*p3^2 - rationalize(BigInt, 0.7692994416588633, tol=tol)*g3^2*p3^2 - rationalize(BigInt, 0.19466145678474384, tol=tol)*q0 - rationalize(BigInt, 0.7153041427190404, tol=tol)*p1*q0 - rationalize(BigInt, 1.3528776260043915, tol=tol)*p2*q0 - rationalize(BigInt, 1.7070452538121381, tol=tol)*p3*q0 - rationalize(BigInt, 1.0516635822669562, tol=tol)*q1 + rationalize(BigInt, 1.2244185478631853, tol=tol)*p1*q1 - rationalize(BigInt, 0.05844567698552443, tol=tol)*p2*q1 - rationalize(BigInt, 0.37706149953585283, tol=tol)*p3*q1 + rationalize(BigInt, 0.580102254517945, tol=tol)*q2 + rationalize(BigInt, 1.2898860704586343, tol=tol)*p1*q2 - rationalize(BigInt, 0.6655948497180294, tol=tol)*p2*q2 + rationalize(BigInt, 0.697758704890495, tol=tol)*p3*q2 - rationalize(BigInt, 0.042921436747585445, tol=tol)*q3 + rationalize(BigInt, 0.5172073855756967, tol=tol)*p1*q3 + rationalize(BigInt, 0.6917094054122289, tol=tol)*p2*q3 - rationalize(BigInt, 1.4579672250860476, tol=tol)*p3*q3,
    rationalize(BigInt, 0.16011034303688113, tol=tol)*g1^2 - rationalize(BigInt, 0.9005468824403076, tol=tol)*g1*g2 - rationalize(BigInt, 0.3519015838689263, tol=tol)*g2^2 + rationalize(BigInt, 0.5202586158306898, tol=tol)*g1*g3 + rationalize(BigInt, 0.908682123022068, tol=tol)*g2*g3 - rationalize(BigInt, 0.4464562170645777, tol=tol)*g3^2 - rationalize(BigInt, 0.13844524415679324, tol=tol)*g1^2*p1 + rationalize(BigInt, 1.5568085644333742, tol=tol)*g1*g2*p1 + rationalize(BigInt, 1.6863862382239232, tol=tol)*g2^2*p1 - rationalize(BigInt, 1.7409458121154344, tol=tol)*g1*g3*p1 - rationalize(BigInt, 0.13872356093602894, tol=tol)*g2*g3*p1 - rationalize(BigInt, 0.5159047084859331, tol=tol)*g3^2*p1 - rationalize(BigInt, 0.2741643484200128, tol=tol)*g1^2*p1^2 - rationalize(BigInt, 0.34212012775550327, tol=tol)*g1*g2*p1^2 - rationalize(BigInt, 0.07542436599114127, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 0.37458987278720324, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 0.4782561996467687, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.3495887144111541, tol=tol)*g3^2*p1^2 + rationalize(BigInt, 0.41377445473869573, tol=tol)*g1^2*p2 - rationalize(BigInt, 1.5789383736211624, tol=tol)*g1*g2*p2 + rationalize(BigInt, 1.268319517294935, tol=tol)*g2^2*p2 + rationalize(BigInt, 0.6163793667190677, tol=tol)*g1*g3*p2 - rationalize(BigInt, 0.43374574206406646, tol=tol)*g2*g3*p2 - rationalize(BigInt, 0.2061458017243186, tol=tol)*g3^2*p2 + rationalize(BigInt, 0.14555549639831628, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 1.1674745895517964, tol=tol)*g1*g2*p1*p2 - rationalize(BigInt, 0.9428064489876502, tol=tol)*g2^2*p1*p2 + rationalize(BigInt, 0.0024916775818734295, tol=tol)*g1*g3*p1*p2 + rationalize(BigInt, 0.5291621555283466, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 0.7972509525893339, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.1807885464109201, tol=tol)*g1^2*p2^2 + rationalize(BigInt, 0.9404541869824675, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.5780030515551372, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 1.0257418447585547, tol=tol)*g1*g3*p2^2 + rationalize(BigInt, 0.09251778173989315, tol=tol)*g2*g3*p2^2 + rationalize(BigInt, 0.39721450514421713, tol=tol)*g3^2*p2^2 + rationalize(BigInt, 0.40272988912109214, tol=tol)*g1^2*p3 - rationalize(BigInt, 0.8272484673958682, tol=tol)*g1*g2*p3 + rationalize(BigInt, 1.057139636924469, tol=tol)*g2^2*p3 - rationalize(BigInt, 0.12353226665002319, tol=tol)*g1*g3*p3 - rationalize(BigInt, 2.5741855761862396, tol=tol)*g2*g3*p3 + rationalize(BigInt, 1.560474007685759, tol=tol)*g3^2*p3 - rationalize(BigInt, 0.6150996832616941, tol=tol)*g1^2*p1*p3 + rationalize(BigInt, 0.09937192239106099, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 0.8226042775491553, tol=tol)*g2^2*p1*p3 + rationalize(BigInt, 0.4732438203631739, tol=tol)*g1*g3*p1*p3 + rationalize(BigInt, 1.6946050580334677, tol=tol)*g2*g3*p1*p3 - rationalize(BigInt, 0.20750459428746135, tol=tol)*g3^2*p1*p3 + rationalize(BigInt, 0.7556396990592089, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 1.412614951501404, tol=tol)*g1*g2*p2*p3 - rationalize(BigInt, 0.09676545515565128, tol=tol)*g2^2*p2*p3 - rationalize(BigInt, 0.9781576342585658, tol=tol)*g1*g3*p2*p3 + rationalize(BigInt, 2.5006617995144724, tol=tol)*g2*g3*p2*p3 - rationalize(BigInt, 0.6588742439035575, tol=tol)*g3^2*p2*p3 + rationalize(BigInt, 0.09337580200909272, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 0.5983340592269643, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.6534274175462785, tol=tol)*g2^2*p3^2 + rationalize(BigInt, 0.6511519719713513, tol=tol)*g1*g3*p3^2 - rationalize(BigInt, 0.5707739813866619, tol=tol)*g2*g3*p3^2 - rationalize(BigInt, 0.7468032195553712, tol=tol)*g3^2*p3^2 + rationalize(BigInt, 0.6382474578966228, tol=tol)*q0 - rationalize(BigInt, 1.032036285581197, tol=tol)*p1*q0 - rationalize(BigInt, 1.4759481703093122, tol=tol)*p2*q0 - rationalize(BigInt, 3.0203435337313205, tol=tol)*p3*q0 - rationalize(BigInt, 0.22812438675350769, tol=tol)*q1 - rationalize(BigInt, 0.2157590670168509, tol=tol)*p1*q1 - rationalize(BigInt, 0.1270558344695696, tol=tol)*p2*q1 - rationalize(BigInt, 0.5148593639524484, tol=tol)*p3*q1 + rationalize(BigInt, 1.3667793800860086, tol=tol)*q2 - rationalize(BigInt, 0.06171123442924746, tol=tol)*p1*q2 - rationalize(BigInt, 0.7314954155886625, tol=tol)*p2*q2 + rationalize(BigInt, 0.7189348075213543, tol=tol)*p3*q2 - rationalize(BigInt, 0.902118536026858, tol=tol)*q3 + rationalize(BigInt, 0.43214823742186254, tol=tol)*p1*q3 + rationalize(BigInt, 0.6677624868260497, tol=tol)*p2*q3 + rationalize(BigInt, 0.5162571144422815, tol=tol)*p3*q3,
    rationalize(BigInt, 0.20816475809219404, tol=tol)*g1^2 - rationalize(BigInt, 0.44624795696445435, tol=tol)*g1*g2 + rationalize(BigInt, 0.1573457781818856, tol=tol)*g2^2 + rationalize(BigInt, 0.2432511536576595, tol=tol)*g1*g3 - rationalize(BigInt, 0.4587424991969163, tol=tol)*g2*g3 - rationalize(BigInt, 0.048748564896809544, tol=tol)*g3^2 - rationalize(BigInt, 0.542424299098038, tol=tol)*g1^2*p1 + rationalize(BigInt, 0.412606879197033, tol=tol)*g1*g2*p1 + rationalize(BigInt, 0.026228493490255755, tol=tol)*g2^2*p1 - rationalize(BigInt, 0.7135454436169615, tol=tol)*g1*g3*p1 + rationalize(BigInt, 0.38566725023570736, tol=tol)*g2*g3*p1 + rationalize(BigInt, 0.06641935708182738, tol=tol)*g3^2*p1 + rationalize(BigInt, 0.1946156279601214, tol=tol)*g1^2*p1^2 + rationalize(BigInt, 1.0726514255671113, tol=tol)*g1*g2*p1^2 - rationalize(BigInt, 0.29746151974577967, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 0.3474150051655493, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 1.5662794253637933, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.10284589178565828, tol=tol)*g3^2*p1^2 + rationalize(BigInt, 0.8655281158446179, tol=tol)*g1^2*p2 - rationalize(BigInt, 1.4227007533612923, tol=tol)*g1*g2*p2  - rationalize(BigInt, 0.3561608986253729, tol=tol)*g2^2*p2 - rationalize(BigInt, 0.7002053827479838, tol=tol)*g1*g3*p2 + rationalize(BigInt, 0.1451263721376322, tol=tol)*g2*g3*p2 - rationalize(BigInt, 0.4446427929457582, tol=tol)*g3^2*p2  - rationalize(BigInt, 0.35089234105147404, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 1.8637325747105546, tol=tol)*g1*g2*p1*p2 - rationalize(BigInt, 0.7643910878410862, tol=tol)*g2^2*p1*p2 + rationalize(BigInt, 0.7421389633104346, tol=tol)*g1*g3*p1*p2 + rationalize(BigInt, 0.8043890896223826, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 1.1152834288925604, tol=tol)*g3^2*p1*p2 - rationalize(BigInt, 0.05058692105297476, tol=tol)*g1^2*p2^2 - rationalize(BigInt, 0.8545531093164939, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.25045809562785276, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 1.482438556873845, tol=tol)*g1*g3*p2^2 - rationalize(BigInt, 0.2760311985894717, tol=tol)*g2*g3*p2^2 + rationalize(BigInt, 0.30104501668082756, tol=tol)*g3^2*p2^2 + rationalize(BigInt, 0.41615153726461007, tol=tol)*g1^2*p3 - rationalize(BigInt, 1.6031132124173149, tol=tol)*g1*g2*p3 + rationalize(BigInt, 1.1652768530802575, tol=tol)*g2^2*p3 + rationalize(BigInt, 0.1236694347662175, tol=tol)*g1*g3*p3  - rationalize(BigInt, 0.033510271732486586, tol=tol)*g2*g3*p3 + rationalize(BigInt, 0.6625023868605743, tol=tol)*g3^2*p3 - rationalize(BigInt, 0.06941899872446193, tol=tol)*g1^2*p1*p3 - rationalize(BigInt, 0.5612725019588681, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 1.4835363108262836, tol=tol)*g2^2*p1*p3 - rationalize(BigInt, 0.8310204341509994, tol=tol)*g1*g3*p1*p3 + rationalize(BigInt, 1.3650887611787323, tol=tol)*g2*g3*p1*p3 - rationalize(BigInt, 1.4141173121018216, tol=tol)*g3^2*p1*p3 - rationalize(BigInt, 0.2915853970368523, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 1.2521117933146961, tol=tol)*g1*g2*p2*p3 + rationalize(BigInt, 0.38706376702247, tol=tol)*g2^2*p2*p3 + rationalize(BigInt, 1.2309129178715645, tol=tol)*g1*g3*p2*p3 + rationalize(BigInt, 2.001338697637118, tol=tol)*g2*g3*p2*p3 - rationalize(BigInt, 0.09547836998561768, tol=tol)*g3^2*p2*p3 - rationalize(BigInt, 0.14402870690714664, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 0.2180983162506176, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.5479196153736324, tol=tol)*g2^2*p3^2 + rationalize(BigInt, 1.1350235517082958, tol=tol)*g1*g3*p3^2 - rationalize(BigInt, 1.2902482267743214, tol=tol)*g2*g3*p3^2 - rationalize(BigInt, 0.40389090846648584, tol=tol)*g3^2*p3^2 - rationalize(BigInt, 0.31676197137727014, tol=tol)*q0 + rationalize(BigInt, 0.44977644852595483, tol=tol)*p1*q0 - rationalize(BigInt, 0.06472442427348668, tol=tol)*p2*q0 - rationalize(BigInt, 2.2439307772054415, tol=tol)*p3*q0 - rationalize(BigInt, 0.5547165223690258, tol=tol)*q1 + rationalize(BigInt, 0.23831878651082344, tol=tol)*p1*q1 + rationalize(BigInt, 0.031977776730485255, tol=tol)*p2*q1 + rationalize(BigInt, 0.16687455406564522, tol=tol)*p3*q1 + rationalize(BigInt, 0.9423377906275198, tol=tol)*q2 + rationalize(BigInt, 1.376589178886685, tol=tol)*p1*q2 + rationalize(BigInt, 0.5306523901876015, tol=tol)*p2*q2 + rationalize(BigInt, 0.4754891181933043, tol=tol)*p3*q2 + rationalize(BigInt, 0.09673230093655334, tol=tol)*q3 + rationalize(BigInt, 0.0892904130224598, tol=tol)*p1*q3 + rationalize(BigInt, 0.943551163213123, tol=tol)*p2*q3 - rationalize(BigInt, 1.2527250130712726, tol=tol)*p3*q3,
     -rationalize(BigInt, 0.04095049824628835, tol=tol)*g1^2 + rationalize(BigInt, 0.043116025511842154, tol=tol)*g1*g2 + rationalize(BigInt, 0.003940499198786224, tol=tol)*g2^2 + rationalize(BigInt, 0.7629770334036455, tol=tol)*g1*g3 - rationalize(BigInt, 0.8492350760146794, tol=tol)*g2*g3 - rationalize(BigInt, 0.2784174783424625, tol=tol)*g3^2 + rationalize(BigInt, 0.4281444292173086, tol=tol)*g1^2*p1 - rationalize(BigInt, 0.8828960936117035, tol=tol)*g1*g2*p1 + rationalize(BigInt, 0.5676682886279524, tol=tol)*g2^2*p1 - rationalize(BigInt, 1.3924527881029736, tol=tol)*g1*g3*p1 - rationalize(BigInt, 0.08671339002537767, tol=tol)*g2*g3*p1 - rationalize(BigInt, 0.7256104095984146, tol=tol)*g3^2*p1 - rationalize(BigInt, 0.5422313181564682, tol=tol)*g1^2*p1^2 + rationalize(BigInt, 0.4871946471731439, tol=tol)*g1*g2*p1^2 + rationalize(BigInt, 0.33193585698170985, tol=tol)*g2^2*p1^2 - rationalize(BigInt, 0.10273772900088107, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 0.6624874115365778, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.21029546117475836, tol=tol)*g3^2*p1^2 + rationalize(BigInt, 1.1356911567255628, tol=tol)*g1^2*p2 - rationalize(BigInt, 1.222100685178249, tol=tol)*g1*g2*p2 - rationalize(BigInt, 0.09334002143332033, tol=tol)*g2^2*p2 - rationalize(BigInt, 1.6524959396527132, tol=tol)*g1*g3*p2 + rationalize(BigInt, 1.5569725124184146, tol=tol)*g2*g3*p2 + rationalize(BigInt, 0.7669386068453008, tol=tol)*g3^2*p2 - rationalize(BigInt, 1.5587675264538823, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 1.1317581527003464, tol=tol)*g1*g2*p1*p2 + rationalize(BigInt, 0.01829080736739283, tol=tol)*g2^2*p1*p2 + rationalize(BigInt, 1.0020439814840232, tol=tol)*g1*g3*p1*p2 + rationalize(BigInt, 0.7984049760283556, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 1.5404767190864894, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.4331904808414006, tol=tol)*g1^2*p2^2 - rationalize(BigInt, 0.1389285705830233, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.012183710127155694, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 0.23769385136664706, tol=tol)*g1*g3*p2^2 + rationalize(BigInt, 0.24382837068218804, tol=tol)*g2*g3*p2^2 - rationalize(BigInt, 0.4210067707142449, tol=tol)*g3^2*p2^2 + rationalize(BigInt, 0.21198555697618326, tol=tol)*g1^2*p3 - rationalize(BigInt, 0.01647709154129892, tol=tol)*g1*g2*p3 - rationalize(BigInt, 0.18945218115272705, tol=tol)*g2^2*p3 + rationalize(BigInt, 0.77528301920843, tol=tol)*g1*g3*p3 - rationalize(BigInt, 2.0514046696465, tol=tol)*g2*g3*p3 + rationalize(BigInt, 0.10432028494512646, tol=tol)*g3^2*p3 - rationalize(BigInt, 0.1676311045685404, tol=tol)*g1^2*p1*p3 - rationalize(BigInt, 1.325546381572095, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 0.8078478514339609, tol=tol)*g2^2*p1*p3 - rationalize(BigInt, 0.9669176278885212, tol=tol)*g1*g3*p1*p3 - rationalize(BigInt, 0.36872926793739896, tol=tol)*g2*g3*p1*p3 - rationalize(BigInt, 0.6402167468654205, tol=tol)*g3^2*p1*p3 + rationalize(BigInt, 0.6297207100844667, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 2.067771321161895, tol=tol)*g1*g2*p2*p3  - rationalize(BigInt, 0.17862819697751522, tol=tol)*g2^2*p2*p3 + rationalize(BigInt, 0.020391323549034297, tol=tol)*g1*g3*p2*p3  + rationalize(BigInt, 2.777563669744398, tol=tol)*g2*g3*p2*p3 - rationalize(BigInt, 0.45109251310695153, tol=tol)*g3^2*p2*p3 + rationalize(BigInt, 0.10904083731506761, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 0.3482660765901206, tol=tol)*g1*g2*p3^2 - rationalize(BigInt, 0.31975214685455416, tol=tol)*g2^2*p3^2 + rationalize(BigInt, 0.34043158036752813, tol=tol)*g1*g3*p3^2 - rationalize(BigInt, 0.9063157822187657, tol=tol)*g2*g3*p3^2 + rationalize(BigInt, 0.21071130953948652, tol=tol)*g3^2*p3^2 + rationalize(BigInt, 0.31542747738996463, tol=tol)*q0 - rationalize(BigInt, 0.27020230824684643, tol=tol)*p1*q0 - rationalize(BigInt, 1.8092897421375431, tol=tol)*p2*q0 - rationalize(BigInt, 0.12685366076858268, tol=tol)*p3*q0 - rationalize(BigInt, 0.9721054060313574, tol=tol)*q1 + rationalize(BigInt, 1.4332583965298273, tol=tol)*p1*q1 - rationalize(BigInt, 0.3658292969614953, tol=tol)*p2*q1 - rationalize(BigInt, 0.458292808629767, tol=tol)*p3*q1 + rationalize(BigInt, 1.1057480001700448, tol=tol)*q2 + rationalize(BigInt, 0.649216154064302, tol=tol)*p1*q2 + rationalize(BigInt, 0.1435470147844548, tol=tol)*p2*q2 + rationalize(BigInt, 1.8049686045262234, tol=tol)*p3*q2 + rationalize(BigInt, 0.3619641675513017, tol=tol)*q3 + rationalize(BigInt, 1.0386298649000567, tol=tol)*p1*q3 - rationalize(BigInt, 0.2739870731830222, tol=tol)*p2*q3 - rationalize(BigInt, 0.38992289294835114, tol=tol)*p3*q3,
     -rationalize(BigInt, 0.41615764608945516, tol=tol)*g1^2 - rationalize(BigInt, 1.2331171001793817, tol=tol)*g1*g2 + rationalize(BigInt, 0.10423594498637195, tol=tol)*g2^2 + rationalize(BigInt, 0.4451741240918564, tol=tol)*g1*g3 - rationalize(BigInt, 0.0807794759847403, tol=tol)*g2*g3 + rationalize(BigInt, 0.015584822151867354, tol=tol)*g3^2 + rationalize(BigInt, 0.5169791211840113, tol=tol)*g1^2*p1 - rationalize(BigInt, 0.3281633186673521, tol=tol)*g1*g2*p1 + rationalize(BigInt, 0.10768082059655043, tol=tol)*g2^2*p1 - rationalize(BigInt, 1.78387184821123, tol=tol)*g1*g3*p1 + rationalize(BigInt, 0.1962385955438586, tol=tol)*g2*g3*p1 - rationalize(BigInt, 0.0932755727182936, tol=tol)*g3^2*p1 - rationalize(BigInt, 0.12307658314371513, tol=tol)*g1^2*p1^2 + rationalize(BigInt, 1.3649915585405705, tol=tol)*g1*g2*p1^2 - rationalize(BigInt, 0.01641144275933561, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 1.282522294958988, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 0.13118389677242223, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.13948802590305073, tol=tol)*g3^2*p1^2 - rationalize(BigInt, 0.4784260518169776, tol=tol)*g1^2*p2 - rationalize(BigInt, 2.092134198423298, tol=tol)*g1*g2*p2 - rationalize(BigInt, 0.2652478875380973, tol=tol)*g2^2*p2 - rationalize(BigInt, 0.9322070346912057, tol=tol)*g1*g3*p2 + rationalize(BigInt, 0.28229645793462466, tol=tol)*g2*g3*p2 - rationalize(BigInt, 0.07438003692790207, tol=tol)*g3^2*p2 + rationalize(BigInt, 0.1376969425780227, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 1.3462021315216954, tol=tol)*g1*g2*p1*p2 - rationalize(BigInt, 0.35831006800801096, tol=tol)*g2^2*p1*p2 - rationalize(BigInt, 0.3236078908735904, tol=tol)*g1*g3*p1*p2 - rationalize(BigInt, 0.23312973570099904, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 0.2206131254299883, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.03212463900566726, tol=tol)*g1^2*p2^2 - rationalize(BigInt, 0.543299069419884, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.1081313237618179, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 1.8348908578280814, tol=tol)*g1*g3*p2^2 - rationalize(BigInt, 0.3379484782876818, tol=tol)*g2*g3*p2^2 + rationalize(BigInt, 0.07600668475615065, tol=tol)*g3^2*p2^2 + rationalize(BigInt, 0.6002572265406737, tol=tol)*g1^2*p3 + rationalize(BigInt, 0.1802055521069689, tol=tol)*g1*g2*p3 + rationalize(BigInt, 1.5644331744196656, tol=tol)*g2^2*p3 - rationalize(BigInt, 0.5267676534246675, tol=tol)*g1*g3*p3 - rationalize(BigInt, 1.676157338774887, tol=tol)*g2*g3*p3 + rationalize(BigInt, 0.4008307913692461, tol=tol)*g3^2*p3 - rationalize(BigInt, 0.5875157312206984, tol=tol)*g1^2*p1*p3 - rationalize(BigInt, 0.26794777918572443, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 1.8161766879761405, tol=tol)*g2^2*p1*p3 + rationalize(BigInt, 0.08217873131336825, tol=tol)*g1*g3*p1*p3 + rationalize(BigInt, 0.41625228622759664, tol=tol)*g2*g3*p1*p3 - rationalize(BigInt, 1.228660956755442, tol=tol)*g3^2*p1*p3 + rationalize(BigInt, 0.801736157469905, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 0.06119905069237429, tol=tol)*g1*g2*p2*p3 + rationalize(BigInt, 0.5320140756032581, tol=tol)*g2^2*p2*p3 - rationalize(BigInt, 0.5120321678148483, tol=tol)*g1*g3*p2*p3 + rationalize(BigInt, 1.3994959245799465, tol=tol)*g2*g3*p2*p3 - rationalize(BigInt, 1.3337502330731632, tol=tol)*g3^2*p2*p3 + rationalize(BigInt, 0.09095194413804787, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 0.8216924891206866, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.12454276652115351, tol=tol)*g2^2*p3^2 + rationalize(BigInt, 0.5523685628690934, tol=tol)*g1*g3*p3^2 + rationalize(BigInt, 0.20676458151525956, tol=tol)*g2*g3*p3^2 - rationalize(BigInt, 0.2154947106592014, tol=tol)*g3^2*p3^2 + rationalize(BigInt, 0.29633687895121585, tol=tol)*q0 - rationalize(BigInt, 0.5313843690622682, tol=tol)*p1*q0 + rationalize(BigInt, 0.8180539762829769, tol=tol)*p2*q0 - rationalize(BigInt, 2.5655211923295855, tol=tol)*p3*q0 - rationalize(BigInt, 1.2602089223582702, tol=tol)*q1 + rationalize(BigInt, 1.1609512634985952, tol=tol)*p1*q1 - rationalize(BigInt, 1.5289512288758575, tol=tol)*p2*q1 - rationalize(BigInt, 0.17424561505966216, tol=tol)*p3*q1 + rationalize(BigInt, 0.10364901876603111, tol=tol)*q2 - rationalize(BigInt, 0.013973408764994696, tol=tol)*p1*q2 - rationalize(BigInt, 0.30080374272031296, tol=tol)*p2*q2 + rationalize(BigInt, 1.5473286276087392, tol=tol)*p3*q2 - rationalize(BigInt, 0.04276145686338927, tol=tol)*q3 + rationalize(BigInt, 0.12497482938060817, tol=tol)*p1*q3 + rationalize(BigInt, 0.06337974500071619, tol=tol)*p2*q3 - rationalize(BigInt, 1.1223229062458282, tol=tol)*p3*q3,
    rationalize(BigInt, 0.0995239208560676, tol=tol)*g1^2 - rationalize(BigInt, 0.4069835646127759, tol=tol)*g1*g2 + rationalize(BigInt, 0.2822661867004282, tol=tol)*g2^2 + rationalize(BigInt, 0.8451107786883508, tol=tol)*g1*g3 - rationalize(BigInt, 0.5316708978452792, tol=tol)*g2*g3 - rationalize(BigInt, 0.8798108963881374, tol=tol)*g3^2 - rationalize(BigInt, 0.08841286657967666, tol=tol)*g1^2*p1 - rationalize(BigInt, 0.5437440767827946, tol=tol)*g1*g2*p1 + rationalize(BigInt, 0.8310187547943032, tol=tol)*g2^2*p1 - rationalize(BigInt, 2.300828301805621, tol=tol)*g1*g3*p1 + rationalize(BigInt, 1.7576212044508612, tol=tol)*g2*g3*p1 + rationalize(BigInt, 0.8123567734335738, tol=tol)*g3^2*p1 - rationalize(BigInt, 0.06820768665581543, tol=tol)*g1^2*p1^2 + rationalize(BigInt, 1.3043060864108098, tol=tol)*g1*g2*p1^2 - rationalize(BigInt, 0.009054065062096384, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 1.3077574118576432, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 0.06768474224643237, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.07726175171791182, tol=tol)*g3^2*p1^2 + rationalize(BigInt, 0.21964796928975452, tol=tol)*g1^2*p2 - rationalize(BigInt, 0.7213869971382741, tol=tol)*g1*g2*p2 + rationalize(BigInt, 0.06428880062023586, tol=tol)*g2^2*p2 + rationalize(BigInt, 0.6769693907722261, tol=tol)*g1*g3*p2 + rationalize(BigInt, 1.0938889186830445, tol=tol)*g2*g3*p2 + rationalize(BigInt, 1.033778176369894, tol=tol)*g3^2*p2 + rationalize(BigInt, 0.04579106142198317, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 1.6136177673887118, tol=tol)*g1*g2*p1*p2 - rationalize(BigInt, 0.08313705215573315, tol=tol)*g2^2*p1*p2 - rationalize(BigInt, 0.3180870028357637, tol=tol)*g1*g3*p1*p2 - rationalize(BigInt, 0.1323253404673372, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 0.03734599073374998, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.06270407236710761, tol=tol)*g1^2*p2^2 - rationalize(BigInt, 0.19162792586107563, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.009068232422623026, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 1.9648502282295741, tol=tol)*g1*g3*p2^2 - rationalize(BigInt, 0.09696983966017489, tol=tol)*g2*g3*p2^2 - rationalize(BigInt, 0.053635839944484585, tol=tol)*g3^2*p2^2 + rationalize(BigInt, 0.09136896607569936, tol=tol)*g1^2*p3 - rationalize(BigInt, 1.5831271499881143, tol=tol)*g1*g2*p3  + rationalize(BigInt, 1.1586203238882382, tol=tol)*g2^2*p3 + rationalize(BigInt, 1.561944037084699, tol=tol)*g1*g3*p3 + rationalize(BigInt, 0.5373765622781439, tol=tol)*g2*g3*p3 - rationalize(BigInt, 0.7510078085771675, tol=tol)*g3^2*p3 + rationalize(BigInt, 0.03313601431114031, tol=tol)*g1^2*p1*p3 + rationalize(BigInt, 0.999575794109116, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 0.0005673645755714211, tol=tol)*g2^2*p1*p3 - rationalize(BigInt, 1.8768450822200746, tol=tol)*g1*g3*p1*p3 + rationalize(BigInt, 0.15574381848783894, tol=tol)*g2*g3*p1*p3 - rationalize(BigInt, 0.03370337888671173, tol=tol)*g3^2*p1*p3 + rationalize(BigInt, 0.040671941002624216, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 2.5271478321380254, tol=tol)*g1*g2*p2*p3 - rationalize(BigInt, 0.11221443152128503, tol=tol)*g2^2*p2*p3 + rationalize(BigInt, 0.5426249873354888, tol=tol)*g1*g3*p2*p3 + rationalize(BigInt, 0.01844664972668835, tol=tol)*g2*g3*p2*p3 + rationalize(BigInt, 0.07154249051866082, tol=tol)*g3^2*p2*p3 + rationalize(BigInt, 0.005503614288707821, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 1.1126781605497342, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.01812229748471941, tol=tol)*g2^2*p3^2 + rationalize(BigInt, 0.657092816371931, tol=tol)*g1*g3*p3^2 + rationalize(BigInt, 0.029285097413742515, tol=tol)*g2*g3*p3^2 - rationalize(BigInt, 0.02362591177342723, tol=tol)*g3^2*p3^2 + rationalize(BigInt, 0.49802078883164164, tol=tol)*q0 - rationalize(BigInt, 1.5549626616482004, tol=tol)*p1*q0 - rationalize(BigInt, 1.3177149462798845, tol=tol)*p2*q0 - rationalize(BigInt, 0.49898148138677, tol=tol)*p3*q0 - rationalize(BigInt, 0.12338824573791841, tol=tol)*q1 - rationalize(BigInt, 0.06111510706623446, tol=tol)*p1*q1 - rationalize(BigInt, 0.04157014428073322, tol=tol)*p2*q1 - rationalize(BigInt, 0.007996832327295606, tol=tol)*p3*q1 + rationalize(BigInt, 0.3953545453627241, tol=tol)*q2 + rationalize(BigInt, 1.168253359907005, tol=tol)*p1*q2 + rationalize(BigInt, 0.1289788318383644, tol=tol)*p2*q2 + rationalize(BigInt, 1.6166074402033737, tol=tol)*p3*q2 - rationalize(BigInt, 1.1634376524391858, tol=tol)*q3 + rationalize(BigInt, 1.1753697910765462, tol=tol)*p1*q3 + rationalize(BigInt, 1.303746265623441, tol=tol)*p2*q3 - rationalize(BigInt, 0.9550529463247988, tol=tol)*p3*q3,
     -rationalize(BigInt, 0.43946259392041137, tol=tol)*g1^2 - rationalize(BigInt, 0.2794194400312886, tol=tol)*g1*g2 + rationalize(BigInt, 0.28633513817241923, tol=tol)*g2^2 + rationalize(BigInt, 0.4837842932289669, tol=tol)*g1*g3 - rationalize(BigInt, 0.2789447622513483, tol=tol)*g2*g3 + rationalize(BigInt, 0.00840391192868436, tol=tol)*g3^2 - rationalize(BigInt, 1.5744202005495247, tol=tol)*g1^2*p1 - rationalize(BigInt, 0.2861663506179279, tol=tol)*g1*g2*p1 + rationalize(BigInt, 0.7386941869051702, tol=tol)*g2^2*p1 + rationalize(BigInt, 0.5331769017373268, tol=tol)*g1*g3*p1 + rationalize(BigInt, 0.3754313238644043, tol=tol)*g2*g3*p1 - rationalize(BigInt, 0.49557272342979514, tol=tol)*g3^2*p1 - rationalize(BigInt, 0.1697121445535322, tol=tol)*g1^2*p1^2 + rationalize(BigInt, 0.17649121044520624, tol=tol)*g1*g2*p1^2 + rationalize(BigInt, 0.230772384999563, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 1.5487222154552227, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 1.1549190062177301, tol=tol)*g2*g3*p1^2 - rationalize(BigInt, 0.061060240446030825, tol=tol)*g3^2*p1^2 + rationalize(BigInt, 0.3428887153466959, tol=tol)*g1^2*p2 - rationalize(BigInt, 2.0520172954900926, tol=tol)*g1*g2*p2 + rationalize(BigInt, 0.2548434377395737, tol=tol)*g2^2*p2 - rationalize(BigInt, 1.1661156986017325, tol=tol)*g1*g3*p2 + rationalize(BigInt, 1.1418016311188108, tol=tol)*g2*g3*p2 - rationalize(BigInt, 0.03717326023724154, tol=tol)*g3^2*p2 - rationalize(BigInt, 1.1298456541690676, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 2.0395244570443807, tol=tol)*g1*g2*p1*p2 - rationalize(BigInt, 0.46158802071478294, tol=tol)*g2^2*p1*p2 - rationalize(BigInt, 0.7937227851031279, tol=tol)*g1*g3*p1*p2 + rationalize(BigInt, 1.3925863980743391, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 1.5914336748838505, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.5461657368466324, tol=tol)*g1^2*p2^2 - rationalize(BigInt, 0.40076456254394, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.5789448963236655, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 1.2606913391501393, tol=tol)*g1*g3*p2^2 - rationalize(BigInt, 0.8839870841626388, tol=tol)*g2*g3*p2^2 + rationalize(BigInt, 0.03277915947703309, tol=tol)*g3^2*p2^2 - rationalize(BigInt, 0.8186969313136523, tol=tol)*g1^2*p3 - rationalize(BigInt, 0.02925076747851844, tol=tol)*g1*g2*p3 + rationalize(BigInt, 0.7392135200738452, tol=tol)*g2^2*p3 + rationalize(BigInt, 0.12848993732517067, tol=tol)*g1*g3*p3 - rationalize(BigInt, 0.9246721961161062, tol=tol)*g2*g3*p3 + rationalize(BigInt, 0.34711806557393554, tol=tol)*g3^2*p3 - rationalize(BigInt, 1.3114980542204153, tol=tol)*g1^2*p1*p3 + rationalize(BigInt, 0.7458693987541652, tol=tol)*g1*g2*p1*p3  + rationalize(BigInt, 1.309443056027262, tol=tol)*g2^2*p1*p3 - rationalize(BigInt, 0.5609050703954811, tol=tol)*g1*g3*p1*p3 - rationalize(BigInt, 0.560863866588139, tol=tol)*g2*g3*p1*p3 + rationalize(BigInt, 0.0020549981931533984, tol=tol)*g3^2*p1*p3 + rationalize(BigInt, 0.21040044489655793, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 2.094674294769259, tol=tol)*g1*g2*p2*p3 + rationalize(BigInt, 0.9014935747766074, tol=tol)*g2^2*p2*p3 - rationalize(BigInt, 0.7741828342315165, tol=tol)*g1*g3*p2*p3 + rationalize(BigInt, 1.0179705141740856, tol=tol)*g2*g3*p2*p3 - rationalize(BigInt, 1.1118940196731655, tol=tol)*g3^2*p2*p3 - rationalize(BigInt, 0.37645359229310016, tol=tol)*g1^2*p3^2 + rationalize(BigInt, 0.22427335209873375, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.34817251132410243, tol=tol)*g2^2*p3^2 - rationalize(BigInt, 0.2880308763050832, tol=tol)*g1*g3*p3^2 - rationalize(BigInt, 0.27093192205509137, tol=tol)*g2*g3*p3^2 + rationalize(BigInt, 0.028281080968997733, tol=tol)*g3^2*p3^2 + rationalize(BigInt, 0.14472354381930777, tol=tol)*q0 + rationalize(BigInt, 1.3312987370741496, tol=tol)*p1*q0 - rationalize(BigInt, 0.560558892849028, tol=tol)*p2*q0 - rationalize(BigInt, 0.26763465433412853, tol=tol)*p3*q0 + rationalize(BigInt, 0.34474968158768826, tol=tol)*q1 + rationalize(BigInt, 1.1967487089081013, tol=tol)*p1*q1 - rationalize(BigInt, 0.5416788729536276, tol=tol)*p2*q1 + rationalize(BigInt, 0.35732161373948407, tol=tol)*p3*q1 + rationalize(BigInt, 0.4086834090152563, tol=tol)*q2 + rationalize(BigInt, 0.9057045158335053, tol=tol)*p1*q2 + rationalize(BigInt, 0.790785041402257, tol=tol)*p2*q2 + rationalize(BigInt, 0.25330644272284614, tol=tol)*p3*q2 - rationalize(BigInt, 0.3854162586341408, tol=tol)*q3 - rationalize(BigInt, 0.04738889206888738, tol=tol)*p1*q3 + rationalize(BigInt, 1.2360880142692436, tol=tol)*p2*q3 - rationalize(BigInt, 0.03146357181747231, tol=tol)*p3*q3,
     -rationalize(BigInt, 0.3364985136694329, tol=tol)*g1^2 - rationalize(BigInt, 0.6220713377668979, tol=tol)*g1*g2  + rationalize(BigInt, 0.30253663989085705, tol=tol)*g2^2 + rationalize(BigInt, 0.09691241616903694, tol=tol)*g1*g3 - rationalize(BigInt, 0.23475287081412524, tol=tol)*g2*g3 + rationalize(BigInt, 0.03759204880218851, tol=tol)*g3^2 - rationalize(BigInt, 1.0655428995245861, tol=tol)*g1^2*p1 + rationalize(BigInt, 1.2015525036653503, tol=tol)*g1*g2*p1 + rationalize(BigInt, 0.7601103825805371, tol=tol)*g2^2*p1 - rationalize(BigInt, 0.19738187400246698, tol=tol)*g1*g3*p1 + rationalize(BigInt, 0.6151569823342218, tol=tol)*g2*g3*p1 - rationalize(BigInt, 0.19896444557080828, tol=tol)*g3^2*p1 + rationalize(BigInt, 0.12867336914719074, tol=tol)*g1^2*p1^2 + rationalize(BigInt, 1.2422164684570598, tol=tol)*g1*g2*p1^2 - rationalize(BigInt, 0.32698199621342966, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 0.988903111633841, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 0.015998020944011523, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.19830862706623895, tol=tol)*g3^2*p1^2 - rationalize(BigInt, 0.36506543203523656, tol=tol)*g1^2*p2 - rationalize(BigInt, 1.8191248166962, tol=tol)*g1*g2*p2 - rationalize(BigInt, 0.04702353062492582, tol=tol)*g2^2*p2 + rationalize(BigInt, 0.22650992191881464, tol=tol)*g1*g3*p2 + rationalize(BigInt, 0.6757394159618821, tol=tol)*g2*g3*p2 - rationalize(BigInt, 0.20784006514510822, tol=tol)*g3^2*p2 - rationalize(BigInt, 1.3823450011564475, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 0.10921823662847378, tol=tol)*g1*g2*p1*p2 + rationalize(BigInt, 0.9411204746791397, tol=tol)*g2^2*p1*p2 + rationalize(BigInt, 0.6416176172245012, tol=tol)*g1*g3*p1*p2 + rationalize(BigInt, 0.6050178988297725, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 0.4412245264773078, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.06735301456610396, tol=tol)*g1^2*p2^2 - rationalize(BigInt, 0.8585408362609629, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.3090762258515019, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 0.4001418483728859, tol=tol)*g1*g3*p2^2 + rationalize(BigInt, 0.4779675484647562, tol=tol)*g2*g3*p2^2 + rationalize(BigInt, 0.24172321128539795, tol=tol)*g3^2*p2^2 - rationalize(BigInt, 0.5137007261636957, tol=tol)*g1^2*p3 - rationalize(BigInt, 0.9820327320175907, tol=tol)*g1*g2*p3 + rationalize(BigInt, 0.926095293449879, tol=tol)*g2^2*p3 - rationalize(BigInt, 0.7003487284409724, tol=tol)*g1*g3*p3 - rationalize(BigInt, 0.9430374895060032, tol=tol)*g2*g3*p3 + rationalize(BigInt, 0.08868217302851986, tol=tol)*g3^2*p3 - rationalize(BigInt, 0.8028419556001375, tol=tol)*g1^2*p1*p3 + rationalize(BigInt, 1.4979362188623169, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 0.6798639272820932, tol=tol)*g2^2*p1*p3 - rationalize(BigInt, 1.0625669164194511, tol=tol)*g1*g3*p1*p3  + rationalize(BigInt, 1.6499737530190743, tol=tol)*g2*g3*p1*p3 + rationalize(BigInt, 0.12297802831804432, tol=tol)*g3^2*p1*p3 - rationalize(BigInt, 0.2829951083769471, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 1.67180387717421, tol=tol)*g1*g2*p2*p3 + rationalize(BigInt, 0.22862840829725906, tol=tol)*g2^2*p2*p3 - rationalize(BigInt, 0.45573353517584536, tol=tol)*g1*g3*p2*p3 - rationalize(BigInt, 0.030399560249614144, tol=tol)*g2*g3*p2*p3 + rationalize(BigInt, 0.05436670007968804, tol=tol)*g3^2*p2*p3 - rationalize(BigInt, 0.1960263837132947, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 0.3836756321960968, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.6360582220649316, tol=tol)*g2^2*p3^2 - rationalize(BigInt, 0.5887612632609552, tol=tol)*g1*g3*p3^2 - rationalize(BigInt, 0.49396556940876774, tol=tol)*g2*g3*p3^2 - rationalize(BigInt, 0.4400318383516369, tol=tol)*g3^2*p3^2 - rationalize(BigInt, 0.0036301750236126602, tol=tol)*q0 + rationalize(BigInt, 0.5043969625148573, tol=tol)*p1*q0 + rationalize(BigInt, 0.6199290278052706, tol=tol)*p2*q0 - rationalize(BigInt, 0.5010767403147031, tol=tol)*p3*q0 + rationalize(BigInt, 0.27845755063810307, tol=tol)*q1 + rationalize(BigInt, 0.9141858220671707, tol=tol)*p1*q1 - rationalize(BigInt, 0.04474642022970836, tol=tol)*p2*q1 + rationalize(BigInt, 0.21510107355697375, tol=tol)*p3*q1 + rationalize(BigInt, 0.6261152163591358, tol=tol)*q2 - rationalize(BigInt, 0.23441023070199246, tol=tol)*p1*q2  + rationalize(BigInt, 0.5860551333365724, tol=tol)*p2*q2 + rationalize(BigInt, 0.6514610848268766, tol=tol)*p3*q2 - rationalize(BigInt, 0.1414398428437503, tol=tol)*q3 + rationalize(BigInt, 0.18837597417330545, tol=tol)*p1*q3 + rationalize(BigInt, 0.23529746996939166, tol=tol)*p2*q3 + rationalize(BigInt, 0.34502969365002717, tol=tol)*p3*q3,
     -rationalize(BigInt, 0.4997153826800627, tol=tol)*g1^2 + rationalize(BigInt, 1.1379854833548109, tol=tol)*g1*g2 - rationalize(BigInt, 0.6474309248194395, tol=tol)*g2^2 + rationalize(BigInt, 0.5441177673182162, tol=tol)*g1*g3  - rationalize(BigInt, 0.6356209362624222, tol=tol)*g2*g3 - rationalize(BigInt, 0.00249167421336129, tol=tol)*g3^2 - rationalize(BigInt, 1.5544041627556424, tol=tol)*g1^2*p1 + rationalize(BigInt, 2.3739950333676623, tol=tol)*g1*g2*p1  - rationalize(BigInt, 0.626498203472872, tol=tol)*g2^2*p1 - rationalize(BigInt, 1.251167893581656, tol=tol)*g1*g3*p1 + rationalize(BigInt, 1.014528145879352, tol=tol)*g2*g3*p1 - rationalize(BigInt, 0.07360314669929015, tol=tol)*g3^2*p1 - rationalize(BigInt, 0.11700211538167525, tol=tol)*g1^2*p1^2 + rationalize(BigInt, 1.4292168664824174, tol=tol)*g1*g2*p1^2 - rationalize(BigInt, 0.08561887506760898, tol=tol)*g2^2*p1^2 + rationalize(BigInt, 0.11756149787203274, tol=tol)*g1*g3*p1^2 + rationalize(BigInt, 1.2753331275027957, tol=tol)*g2*g3*p1^2 + rationalize(BigInt, 0.20262099044928422, tol=tol)*g3^2*p1^2 - rationalize(BigInt, 0.5322121916934887, tol=tol)*g1^2*p2 - rationalize(BigInt, 0.5705836762208184, tol=tol)*g1*g2*p2 + rationalize(BigInt, 1.338552997747994, tol=tol)*g2^2*p2 + rationalize(BigInt, 0.1626686039427204, tol=tol)*g1*g3*p2 + rationalize(BigInt, 0.8255938648528275, tol=tol)*g2*g3*p2 - rationalize(BigInt, 0.2877403739032339, tol=tol)*g3^2*p2 - rationalize(BigInt, 1.7484788356484942, tol=tol)*g1^2*p1*p2 - rationalize(BigInt, 0.7462669895235802, tol=tol)*g1*g2*p1*p2 + rationalize(BigInt, 0.9024930142580471, tol=tol)*g2^2*p1*p2 - rationalize(BigInt, 0.836428661723292, tol=tol)*g1*g3*p1*p2 + rationalize(BigInt, 0.10315211869079574, tol=tol)*g2*g3*p1*p2 + rationalize(BigInt, 0.8459858213904471, tol=tol)*g3^2*p1*p2 + rationalize(BigInt, 0.0524174914208927, tol=tol)*g1^2*p2^2 - rationalize(BigInt, 1.0349080883999526, tol=tol)*g1*g2*p2^2 - rationalize(BigInt, 0.4456138068634661, tol=tol)*g2^2*p2^2 - rationalize(BigInt, 0.8702119648626987, tol=tol)*g1*g3*p2^2 + rationalize(BigInt, 0.1347633798985455, tol=tol)*g2*g3*p2^2 + rationalize(BigInt, 0.3931963154425734, tol=tol)*g3^2*p2^2 + rationalize(BigInt, 0.7263228121576023, tol=tol)*g1^2*p3 - rationalize(BigInt, 2.2210658961930227, tol=tol)*g1*g2*p3 + rationalize(BigInt, 1.6108892968930464, tol=tol)*g2^2*p3 - rationalize(BigInt, 1.0504867050871953, tol=tol)*g1*g3*p3 + rationalize(BigInt, 1.0659099359741329, tol=tol)*g2*g3*p3 + rationalize(BigInt, 0.5672833993063003, tol=tol)*g3^2*p3 - rationalize(BigInt, 0.06811413086649659, tol=tol)*g1^2*p1*p3 - rationalize(BigInt, 0.31001429706461026, tol=tol)*g1*g2*p1*p3 + rationalize(BigInt, 1.4153449120254278, tol=tol)*g2^2*p1*p3 - rationalize(BigInt, 1.8526453839326555, tol=tol)*g1*g3*p1*p3 + rationalize(BigInt, 1.3321822393355465, tol=tol)*g2*g3*p1*p3 - rationalize(BigInt, 1.3472307811589312, tol=tol)*g3^2*p1*p3 + rationalize(BigInt, 0.8918357940981352, tol=tol)*g1^2*p2*p3 - rationalize(BigInt, 1.500420306224534, tol=tol)*g1*g2*p2*p3 - rationalize(BigInt, 0.43633002867072346, tol=tol)*g2^2*p2*p3 - rationalize(BigInt, 1.087845994748532, tol=tol)*g1*g3*p2*p3 - rationalize(BigInt, 0.15102803926606495, tol=tol)*g2*g3*p2*p3 - rationalize(BigInt, 0.4555057654274117, tol=tol)*g3^2*p2*p3 + rationalize(BigInt, 0.06458462396078254, tol=tol)*g1^2*p3^2 - rationalize(BigInt, 0.39430877808246473, tol=tol)*g1*g2*p3^2 + rationalize(BigInt, 0.5312326819310751, tol=tol)*g2^2*p3^2 + rationalize(BigInt, 0.752650466990666, tol=tol)*g1*g3*p3^2 - rationalize(BigInt, 1.4100965074013412, tol=tol)*g2*g3*p3^2 - rationalize(BigInt, 0.5958173058918577, tol=tol)*g3^2*p3^2 + rationalize(BigInt, 1.1496379817128635, tol=tol)*q0 + rationalize(BigInt, 2.2545055129278047, tol=tol)*p1*q0 - rationalize(BigInt, 0.5186004321512713, tol=tol)*p2*q0 - rationalize(BigInt, 2.904495508356949, tol=tol)*p3*q0 + rationalize(BigInt, 0.3811371457709, tol=tol)*q1 + rationalize(BigInt, 1.1561384774363503, tol=tol)*p1*q1 - rationalize(BigInt, 0.03459051968276095, tol=tol)*p2*q1 - rationalize(BigInt, 0.5860104174434722, tol=tol)*p3*q1 - rationalize(BigInt, 0.42262305475993184, tol=tol)*q2 - rationalize(BigInt, 0.06960283700014054, tol=tol)*p1*q2 + rationalize(BigInt, 0.6975228406191866, tol=tol)*p2*q2 + rationalize(BigInt, 1.1761445768437886, tol=tol)*p3*q2 - rationalize(BigInt, 0.4132504788764938, tol=tol)*q3 + rationalize(BigInt, 1.0477095469315294, tol=tol)*p1*q3 + rationalize(BigInt, 0.5581768684334514, tol=tol)*p2*q3 + rationalize(BigInt, 0.4360576481532085, tol=tol)*p3*q3,
    -1. + rationalize(BigInt, 0.9336143308746049, tol=tol)*g1 + rationalize(BigInt, 1.1781580271766483, tol=tol)*g2 + rationalize(BigInt, 0.551650235747964, tol=tol)*g3
    ]
    equations
end

# Source:
# https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/master/src/systems.jl
#
# chandran(10): 
# Critical pair degree sequence: 3,2,3,4,5,6,7,8,9,10,11
function chandran(n; tol=0, np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, H = np.polynomial_ring(k, ["H$i" for i in 1:n], internal_ordering=internal_ordering)
    c = rationalize(BigInt, 0.51234, tol=tol)
    c = k(numerator(c)) // k(denominator(c))
    eqs = [(2n*H[i] - c*H[i]*(1 + sum(k(i) // (j+i) * H[j] for j=1:(n-1))) - 2n) for i=1:n]
    eqs
end

# Source:
# https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/master/src/systems.jl
function boon(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (s1, g1, s2, g2, C1, C2) = np.polynomial_ring(k, ["s1", "g1", "s2", "g2", "C1", "C2"], internal_ordering=internal_ordering)
    eqs = [
        s1^2+g1^2 - 1,
        s2^2+g2^2 - 1,
        C1*g1^3+C2*g2^3 - k(12)//10,
        C1*s1^3+C2*s2^3 - k(12)//10,
        C1*g1^2*s1+C2*g2^2*s2 - k(7)//10,
        C1*g1*s1^2+C2*g2*s2^2 - k(7)//10
    ]
    eqs
end

# Source:
# https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/e04087ef08cf91ffafd88546c4d4ccb25613a3c7/src/systems.jl#L151
function ipp(; tol=0, np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    R, (x1, x2, x3, x4, x5, x6, x7, x8) = np.polynomial_ring(np.QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"], internal_ordering=internal_ordering)
    sys = [
        x1^2+x2^2-1,
        x3^2+x4^2-1,
        x5^2+x6^2-1,
        x7^2+x8^2-1,
        (-rationalize(BigInt, 2.4915068E-01, tol=tol)*x1*x3+ rationalize(BigInt, 1.6091354E+00, tol=tol)*x1*x4+ rationalize(BigInt, 2.7942343E-01, tol=tol) *x2*x3+ rationalize(BigInt, 1.4348016E+00, tol=tol)*x2*x4) + (rationalize(BigInt, 4.0026384E-01, tol=tol)*x5*x8-rationalize(BigInt, 8.0052768E-01, tol=tol)*x6*x7+ rationalize(BigInt, 7.4052388E-02, tol=tol)*x1-rationalize(BigInt, 8.3050031E-02, tol=tol)*x2) - (rationalize(BigInt, 3.8615961E-01, tol=tol)*x3-rationalize(BigInt, 7.5526603E-01, tol=tol)*x4+ rationalize(BigInt, 5.0420168E-01, tol=tol)*x5 -rationalize(BigInt, 1.0916287E+00, tol=tol)*x6+ rationalize(BigInt, 4.0026384E-01, tol=tol)*x8) + rationalize(BigInt, 4.920729E-02, tol=tol),
        (rationalize(BigInt, 1.2501635E-01, tol=tol)*x1*x3-rationalize(BigInt, 6.8660736E-01, tol=tol)*x1*x4-rationalize(BigInt, 1.1922812E-01, tol=tol)* x2*x3-rationalize(BigInt, 7.1994047E-01, tol=tol)*x2*x4) - (rationalize(BigInt, 4.3241927E-01, tol=tol)*x5*x7-rationalize(BigInt, 8.6483855E-01, tol=tol)*x6*x8-rationalize(BigInt, 3.715727E-02, tol=tol)*x1+ rationalize(BigInt, 3.5436896E-02, tol=tol)*x2)+ rationalize(BigInt, 8.5383482E-02, tol=tol)*x3-rationalize(BigInt, 3.9251967E-02, tol=tol)*x5-rationalize(BigInt, 4.3241927E-01, tol=tol)*x7+ rationalize(BigInt, 1.387301E-02, tol=tol),
        (-rationalize(BigInt, 6.3555007E-01, tol=tol)*x1*x3-rationalize(BigInt, 1.1571992E-01, tol=tol)*x1*x4-rationalize(BigInt, 6.6640448E-01, tol=tol) *x2*x3) + (rationalize(BigInt, 1.1036211E-01, tol=tol)*x2*x4+ rationalize(BigInt, 2.9070203E-01, tol=tol)*x5*x7+ rationalize(BigInt, 1.2587767E+00, tol=tol)*x5*x8)- (rationalize(BigInt, 6.2938836E-01, tol=tol)*x6*x7+ rationalize(BigInt, 5.8140406E-01, tol=tol)*x6*x8+ rationalize(BigInt, 1.9594662E-01, tol=tol)*x1)- (rationalize(BigInt, 1.2280342E+00, tol=tol)*x2-rationalize(BigInt, 7.9034221E-02, tol=tol)*x4+ rationalize(BigInt, 2.6387877E-02, tol=tol)*x5)- rationalize(BigInt, 5.713143E-02, tol=tol)*x6-rationalize(BigInt, 1.1628081E+00, tol=tol)*x7+rationalize(BigInt, 1.2587767E+00, tol=tol)*x8+ rationalize(BigInt, 2.162575E+00, tol=tol),
        (rationalize(BigInt, 1.4894773E+00, tol=tol)*x1*x3+ rationalize(BigInt, 2.3062341E-01, tol=tol)*x1*x4+ rationalize(BigInt, 1.3281073E+00, tol=tol)*x2*x3)-(rationalize(BigInt, 2.5864503E-01, tol=tol)*x2*x4+ rationalize(BigInt, 1.165172E+00, tol=tol)*x5*x7-rationalize(BigInt, 2.6908494E-01, tol=tol)*x5*x8)+ (rationalize(BigInt, 5.3816987E-01, tol=tol)*x6*x7+ rationalize(BigInt, 5.8258598E-01, tol=tol)*x6*x8-rationalize(BigInt, 2.0816985E-01, tol=tol)*x1)+(rationalize(BigInt, 2.686832E+00, tol=tol)*x2-rationalize(BigInt, 6.9910317E-01, tol=tol)*x3+ rationalize(BigInt, 3.5744413E-01, tol=tol)*x4)+ rationalize(BigInt, 1.2499117E+00, tol=tol)*x5+ rationalize(BigInt, 1.467736E+00, tol=tol)*x6+ rationalize(BigInt, 1.165172E+00, tol=tol)*x7+ rationalize(BigInt, 1.10763397E+00, tol=tol)*x8-rationalize(BigInt, 6.9686809E-01, tol=tol)
    ]
    sys = map(f -> AbstractAlgebra.map_coefficients(c -> k(BigInt(numerator(c))) // k(BigInt(denominator(c))), f), sys)
    sys
end

###

# Source: https://web.archive.org/web/20201202185136/http://www.cecm.sfu.ca/%7Erpearcea/mgb.html
function yang1(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48) = np.polynomial_ring(
        k, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14, :x15, :x16, :x17, :x18, :x19, :x20, :x21, :x22, :x23, :x24, :x25, :x26, :x27, :x28, :x29, :x30, :x31, :x32, :x33, :x34, :x35, :x36, :x37, :x38, :x39, :x40, :x41, :x42, :x43, :x44, :x45, :x46, :x47, :x48], internal_ordering=internal_ordering)
    sys = [
        x21*x45+x22*x46+x23*x47+x24*x48, x17*x45+x18*x46+x19*x47+x20*x48, x13*x45+x14*x46+x15*x47+x16*x48, x10*x46+x11*x47+x12*x48+x45*x9, x45*x5+x46*x6+x47*x7+x48*x8, x1*x45+x2*x46+x3*x47+x4*x48, x21*x41+x22*x42+x23*x43+x24*x44,
        x17*x41+x18*x42+x19*x43+x20*x44, x13*x41+x14*x42+x15*x43+x16*x44, x10*x42+x11*x43+x12*x44+x41*x9, x41*x5+x42*x6+x43*x7+x44*x8, x1*x41+x2*x42+x3*x43+x4*x44, x21*x37+x22*x38+x23*x39+x24*x40, x17*x37+x18*x38+x19*x39+x20*x40,
        x13*x37+x14*x38+x15*x39+x16*x40, x10*x38+x11*x39+x12*x40+x37*x9, x37*x5+x38*x6+x39*x7+x40*x8, x1*x37+x2*x38+x3*x39+x4*x40, x21*x33+x22*x34+x23*x35+x24*x36, x17*x33+x18*x34+x19*x35+x20*x36, x13*x33+x14*x34+x15*x35+x16*x36,
        x10*x34+x11*x35+x12*x36+x33*x9, x33*x5+x34*x6+x35*x7+x36*x8, x1*x33+x2*x34+x3*x35+x36*x4, x21*x29+x22*x30+x23*x31+x24*x32, x17*x29+x18*x30+x19*x31+x20*x32, x13*x29+x14*x30+x15*x31+x16*x32, x10*x30+x11*x31+x12*x32+x29*x9,
        x29*x5+x30*x6+x31*x7+x32*x8, x1*x29+x2*x30+x3*x31+x32*x4, x21*x25+x22*x26+x23*x27+x24*x28, x17*x25+x18*x26+x19*x27+x20*x28, x13*x25+x14*x26+x15*x27+x16*x28, x10*x26+x11*x27+x12*x28+x25*x9, x25*x5+x26*x6+x27*x7+x28*x8, x1*x25+x2*x26+x27*x3+x28*x4,
        x33*x38*x43*x48-x33*x38*x44*x47-x33*x39*x42*x48+x33*x39*x44*x46+x33*x40*x42*x47-x33*x40*x43*x46-x34*x37*x43*x48+x34*x37*x44*x47+x34*x39*x41*x48-x34*x39*x44*x45-x34*x40*x41*x47+x34*x40*x43*x45+x35*x37*x42*x48-x35*x37*x44*x46-x35*x38*x41*x48+x35*x38*x44*x45+x35*x40*x41*x46-x35*x40*x42*x45-x36*x37*x42*x47+x36*x37*x43*x46+x36*x38*x41*x47-x36*x38*x43*x45-x36*x39*x41*x46+x36*x39*x42*x45,
        x29*x38*x43*x48-x29*x38*x44*x47-x29*x39*x42*x48+x29*x39*x44*x46+x29*x40*x42*x47-x29*x40*x43*x46-x30*x37*x43*x48+x30*x37*x44*x47+x30*x39*x41*x48-x30*x39*x44*x45-x30*x40*x41*x47+x30*x40*x43*x45+x31*x37*x42*x48-x31*x37*x44*x46-x31*x38*x41*x48+x31*x38*x44*x45+x31*x40*x41*x46-x31*x40*x42*x45-x32*x37*x42*x47+x32*x37*x43*x46+x32*x38*x41*x47-x32*x38*x43*x45-x32*x39*x41*x46+x32*x39*x42*x45,
        x25*x38*x43*x48-x25*x38*x44*x47-x25*x39*x42*x48+x25*x39*x44*x46+x25*x40*x42*x47-x25*x40*x43*x46-x26*x37*x43*x48+x26*x37*x44*x47+x26*x39*x41*x48-x26*x39*x44*x45-x26*x40*x41*x47+x26*x40*x43*x45+x27*x37*x42*x48-x27*x37*x44*x46-x27*x38*x41*x48+x27*x38*x44*x45+x27*x40*x41*x46-x27*x40*x42*x45-x28*x37*x42*x47+x28*x37*x43*x46+x28*x38*x41*x47-x28*x38*x43*x45-x28*x39*x41*x46+x28*x39*x42*x45,
        x29*x34*x43*x48-x29*x34*x44*x47-x29*x35*x42*x48+x29*x35*x44*x46+x29*x36*x42*x47-x29*x36*x43*x46-x30*x33*x43*x48+x30*x33*x44*x47+x30*x35*x41*x48-x30*x35*x44*x45-x30*x36*x41*x47+x30*x36*x43*x45+x31*x33*x42*x48-x31*x33*x44*x46-x31*x34*x41*x48+x31*x34*x44*x45+x31*x36*x41*x46-x31*x36*x42*x45-x32*x33*x42*x47+x32*x33*x43*x46+x32*x34*x41*x47-x32*x34*x43*x45-x32*x35*x41*x46+x32*x35*x42*x45,
        x25*x34*x43*x48-x25*x34*x44*x47-x25*x35*x42*x48+x25*x35*x44*x46+x25*x36*x42*x47-x25*x36*x43*x46-x26*x33*x43*x48+x26*x33*x44*x47+x26*x35*x41*x48-x26*x35*x44*x45-x26*x36*x41*x47+x26*x36*x43*x45+x27*x33*x42*x48-x27*x33*x44*x46-x27*x34*x41*x48+x27*x34*x44*x45+x27*x36*x41*x46-x27*x36*x42*x45-x28*x33*x42*x47+x28*x33*x43*x46+x28*x34*x41*x47-x28*x34*x43*x45-x28*x35*x41*x46+x28*x35*x42*x45,
        x25*x30*x43*x48-x25*x30*x44*x47-x25*x31*x42*x48+x25*x31*x44*x46+x25*x32*x42*x47-x25*x32*x43*x46-x26*x29*x43*x48+x26*x29*x44*x47+x26*x31*x41*x48-x26*x31*x44*x45-x26*x32*x41*x47+x26*x32*x43*x45+x27*x29*x42*x48-x27*x29*x44*x46-x27*x30*x41*x48+x27*x30*x44*x45+x27*x32*x41*x46-x27*x32*x42*x45-x28*x29*x42*x47+x28*x29*x43*x46+x28*x30*x41*x47-x28*x30*x43*x45-x28*x31*x41*x46+x28*x31*x42*x45,
        x29*x34*x39*x48-x29*x34*x40*x47-x29*x35*x38*x48+x29*x35*x40*x46+x29*x36*x38*x47-x29*x36*x39*x46-x30*x33*x39*x48+x30*x33*x40*x47+x30*x35*x37*x48-x30*x35*x40*x45-x30*x36*x37*x47+x30*x36*x39*x45+x31*x33*x38*x48-x31*x33*x40*x46-x31*x34*x37*x48+x31*x34*x40*x45+x31*x36*x37*x46-x31*x36*x38*x45-x32*x33*x38*x47+x32*x33*x39*x46+x32*x34*x37*x47-x32*x34*x39*x45-x32*x35*x37*x46+x32*x35*x38*x45,
        x25*x34*x39*x48-x25*x34*x40*x47-x25*x35*x38*x48+x25*x35*x40*x46+x25*x36*x38*x47-x25*x36*x39*x46-x26*x33*x39*x48+x26*x33*x40*x47+x26*x35*x37*x48-x26*x35*x40*x45-x26*x36*x37*x47+x26*x36*x39*x45+x27*x33*x38*x48-x27*x33*x40*x46-x27*x34*x37*x48+x27*x34*x40*x45+x27*x36*x37*x46-x27*x36*x38*x45-x28*x33*x38*x47+x28*x33*x39*x46+x28*x34*x37*x47-x28*x34*x39*x45-x28*x35*x37*x46+x28*x35*x38*x45,
        x25*x30*x39*x48-x25*x30*x40*x47-x25*x31*x38*x48+x25*x31*x40*x46+x25*x32*x38*x47-x25*x32*x39*x46-x26*x29*x39*x48+x26*x29*x40*x47+x26*x31*x37*x48-x26*x31*x40*x45-x26*x32*x37*x47+x26*x32*x39*x45+x27*x29*x38*x48-x27*x29*x40*x46-x27*x30*x37*x48+x27*x30*x40*x45+x27*x32*x37*x46-x27*x32*x38*x45-x28*x29*x38*x47+x28*x29*x39*x46+x28*x30*x37*x47-x28*x30*x39*x45-x28*x31*x37*x46+x28*x31*x38*x45,
        x25*x30*x35*x48-x25*x30*x36*x47-x25*x31*x34*x48+x25*x31*x36*x46+x25*x32*x34*x47-x25*x32*x35*x46-x26*x29*x35*x48+x26*x29*x36*x47+x26*x31*x33*x48-x26*x31*x36*x45-x26*x32*x33*x47+x26*x32*x35*x45+x27*x29*x34*x48-x27*x29*x36*x46-x27*x30*x33*x48+x27*x30*x36*x45+x27*x32*x33*x46-x27*x32*x34*x45-x28*x29*x34*x47+x28*x29*x35*x46+x28*x30*x33*x47-x28*x30*x35*x45-x28*x31*x33*x46+x28*x31*x34*x45,
        x29*x34*x39*x44-x29*x34*x40*x43-x29*x35*x38*x44+x29*x35*x40*x42+x29*x36*x38*x43-x29*x36*x39*x42-x30*x33*x39*x44+x30*x33*x40*x43+x30*x35*x37*x44-x30*x35*x40*x41-x30*x36*x37*x43+x30*x36*x39*x41+x31*x33*x38*x44-x31*x33*x40*x42-x31*x34*x37*x44+x31*x34*x40*x41+x31*x36*x37*x42-x31*x36*x38*x41-x32*x33*x38*x43+x32*x33*x39*x42+x32*x34*x37*x43-x32*x34*x39*x41-x32*x35*x37*x42+x32*x35*x38*x41,
        x25*x34*x39*x44-x25*x34*x40*x43-x25*x35*x38*x44+x25*x35*x40*x42+x25*x36*x38*x43-x25*x36*x39*x42-x26*x33*x39*x44+x26*x33*x40*x43+x26*x35*x37*x44-x26*x35*x40*x41-x26*x36*x37*x43+x26*x36*x39*x41+x27*x33*x38*x44-x27*x33*x40*x42-x27*x34*x37*x44+x27*x34*x40*x41+x27*x36*x37*x42-x27*x36*x38*x41-x28*x33*x38*x43+x28*x33*x39*x42+x28*x34*x37*x43-x28*x34*x39*x41-x28*x35*x37*x42+x28*x35*x38*x41,
        x25*x30*x39*x44-x25*x30*x40*x43-x25*x31*x38*x44+x25*x31*x40*x42+x25*x32*x38*x43-x25*x32*x39*x42-x26*x29*x39*x44+x26*x29*x40*x43+x26*x31*x37*x44-x26*x31*x40*x41-x26*x32*x37*x43+x26*x32*x39*x41+x27*x29*x38*x44-x27*x29*x40*x42-x27*x30*x37*x44+x27*x30*x40*x41+x27*x32*x37*x42-x27*x32*x38*x41-x28*x29*x38*x43+x28*x29*x39*x42+x28*x30*x37*x43-x28*x30*x39*x41-x28*x31*x37*x42+x28*x31*x38*x41,
        x25*x30*x35*x44-x25*x30*x36*x43-x25*x31*x34*x44+x25*x31*x36*x42+x25*x32*x34*x43-x25*x32*x35*x42-x26*x29*x35*x44+x26*x29*x36*x43+x26*x31*x33*x44-x26*x31*x36*x41-x26*x32*x33*x43+x26*x32*x35*x41+x27*x29*x34*x44-x27*x29*x36*x42-x27*x30*x33*x44+x27*x30*x36*x41+x27*x32*x33*x42-x27*x32*x34*x41-x28*x29*x34*x43+x28*x29*x35*x42+x28*x30*x33*x43-x28*x30*x35*x41-x28*x31*x33*x42+x28*x31*x34*x41,
        x25*x30*x35*x40-x25*x30*x36*x39-x25*x31*x34*x40+x25*x31*x36*x38+x25*x32*x34*x39-x25*x32*x35*x38-x26*x29*x35*x40+x26*x29*x36*x39+x26*x31*x33*x40-x26*x31*x36*x37-x26*x32*x33*x39+x26*x32*x35*x37+x27*x29*x34*x40-x27*x29*x36*x38-x27*x30*x33*x40+x27*x30*x36*x37+x27*x32*x33*x38-x27*x32*x34*x37-x28*x29*x34*x39+x28*x29*x35*x38+x28*x30*x33*x39-x28*x30*x35*x37-x28*x31*x33*x38+x28*x31*x34*x37,
        -x10*x13*x19*x24+x10*x13*x20*x23+x10*x15*x17*x24-x10*x15*x20*x21-x10*x16*x17*x23+x10*x16*x19*x21+x11*x13*x18*x24-x11*x13*x20*x22-x11*x14*x17*x24+x11*x14*x20*x21+x11*x16*x17*x22-x11*x16*x18*x21-x12*x13*x18*x23+x12*x13*x19*x22+x12*x14*x17*x23-x12*x14*x19*x21-x12*x15*x17*x22+x12*x15*x18*x21+x14*x19*x24*x9-x14*x20*x23*x9-x15*x18*x24*x9+x15*x20*x22*x9+x16*x18*x23*x9-x16*x19*x22*x9,
        -x13*x18*x23*x8+x13*x18*x24*x7+x13*x19*x22*x8-x13*x19*x24*x6-x13*x20*x22*x7+x13*x20*x23*x6+x14*x17*x23*x8-x14*x17*x24*x7-x14*x19*x21*x8+x14*x19*x24*x5+x14*x20*x21*x7-x14*x20*x23*x5-x15*x17*x22*x8+x15*x17*x24*x6+x15*x18*x21*x8-x15*x18*x24*x5-x15*x20*x21*x6+x15*x20*x22*x5+x16*x17*x22*x7-x16*x17*x23*x6-x16*x18*x21*x7+x16*x18*x23*x5+x16*x19*x21*x6-x16*x19*x22*x5,
        x1*x14*x19*x24-x1*x14*x20*x23-x1*x15*x18*x24+x1*x15*x20*x22+x1*x16*x18*x23-x1*x16*x19*x22-x13*x18*x23*x4+x13*x18*x24*x3-x13*x19*x2*x24+x13*x19*x22*x4+x13*x2*x20*x23-x13*x20*x22*x3+x14*x17*x23*x4-x14*x17*x24*x3-x14*x19*x21*x4+x14*x20*x21*x3+x15*x17*x2*x24-x15*x17*x22*x4+x15*x18*x21*x4-x15*x2*x20*x21-x16*x17*x2*x23+x16*x17*x22*x3-x16*x18*x21*x3+x16*x19*x2*x21,
        x10*x17*x23*x8-x10*x17*x24*x7-x10*x19*x21*x8+x10*x19*x24*x5+x10*x20*x21*x7-x10*x20*x23*x5-x11*x17*x22*x8+x11*x17*x24*x6+x11*x18*x21*x8-x11*x18*x24*x5-x11*x20*x21*x6+x11*x20*x22*x5+x12*x17*x22*x7-x12*x17*x23*x6-x12*x18*x21*x7+x12*x18*x23*x5+x12*x19*x21*x6-x12*x19*x22*x5-x18*x23*x8*x9+x18*x24*x7*x9+x19*x22*x8*x9-x19*x24*x6*x9-x20*x22*x7*x9+x20*x23*x6*x9,
        x1*x10*x19*x24-x1*x10*x20*x23-x1*x11*x18*x24+x1*x11*x20*x22+x1*x12*x18*x23-x1*x12*x19*x22+x10*x17*x23*x4-x10*x17*x24*x3-x10*x19*x21*x4+x10*x20*x21*x3+x11*x17*x2*x24-x11*x17*x22*x4+x11*x18*x21*x4-x11*x2*x20*x21-x12*x17*x2*x23+x12*x17*x22*x3-x12*x18*x21*x3+x12*x19*x2*x21-x18*x23*x4*x9+x18*x24*x3*x9-x19*x2*x24*x9+x19*x22*x4*x9+x2*x20*x23*x9-x20*x22*x3*x9,
        x1*x18*x23*x8-x1*x18*x24*x7-x1*x19*x22*x8+x1*x19*x24*x6+x1*x20*x22*x7-x1*x20*x23*x6-x17*x2*x23*x8+x17*x2*x24*x7+x17*x22*x3*x8-x17*x22*x4*x7+x17*x23*x4*x6-x17*x24*x3*x6-x18*x21*x3*x8+x18*x21*x4*x7-x18*x23*x4*x5+x18*x24*x3*x5+x19*x2*x21*x8-x19*x2*x24*x5-x19*x21*x4*x6+x19*x22*x4*x5-x2*x20*x21*x7+x2*x20*x23*x5+x20*x21*x3*x6-x20*x22*x3*x5,
        x10*x13*x23*x8-x10*x13*x24*x7-x10*x15*x21*x8+x10*x15*x24*x5+x10*x16*x21*x7-x10*x16*x23*x5-x11*x13*x22*x8+x11*x13*x24*x6+x11*x14*x21*x8-x11*x14*x24*x5-x11*x16*x21*x6+x11*x16*x22*x5+x12*x13*x22*x7-x12*x13*x23*x6-x12*x14*x21*x7+x12*x14*x23*x5+x12*x15*x21*x6-x12*x15*x22*x5-x14*x23*x8*x9+x14*x24*x7*x9+x15*x22*x8*x9-x15*x24*x6*x9-x16*x22*x7*x9+x16*x23*x6*x9,
        x1*x10*x15*x24-x1*x10*x16*x23-x1*x11*x14*x24+x1*x11*x16*x22+x1*x12*x14*x23-x1*x12*x15*x22+x10*x13*x23*x4-x10*x13*x24*x3-x10*x15*x21*x4+x10*x16*x21*x3+x11*x13*x2*x24-x11*x13*x22*x4+x11*x14*x21*x4-x11*x16*x2*x21-x12*x13*x2*x23+x12*x13*x22*x3-x12*x14*x21*x3+x12*x15*x2*x21-x14*x23*x4*x9+x14*x24*x3*x9-x15*x2*x24*x9+x15*x22*x4*x9+x16*x2*x23*x9-x16*x22*x3*x9,
        x1*x14*x23*x8-x1*x14*x24*x7-x1*x15*x22*x8+x1*x15*x24*x6+x1*x16*x22*x7-x1*x16*x23*x6-x13*x2*x23*x8+x13*x2*x24*x7+x13*x22*x3*x8-x13*x22*x4*x7+x13*x23*x4*x6-x13*x24*x3*x6-x14*x21*x3*x8+x14*x21*x4*x7-x14*x23*x4*x5+x14*x24*x3*x5+x15*x2*x21*x8-x15*x2*x24*x5-x15*x21*x4*x6+x15*x22*x4*x5-x16*x2*x21*x7+x16*x2*x23*x5+x16*x21*x3*x6-x16*x22*x3*x5,
        x1*x10*x23*x8-x1*x10*x24*x7-x1*x11*x22*x8+x1*x11*x24*x6+x1*x12*x22*x7-x1*x12*x23*x6-x10*x21*x3*x8+x10*x21*x4*x7-x10*x23*x4*x5+x10*x24*x3*x5+x11*x2*x21*x8-x11*x2*x24*x5-x11*x21*x4*x6+x11*x22*x4*x5-x12*x2*x21*x7+x12*x2*x23*x5+x12*x21*x3*x6-x12*x22*x3*x5-x2*x23*x8*x9+x2*x24*x7*x9+x22*x3*x8*x9-x22*x4*x7*x9+x23*x4*x6*x9-x24*x3*x6*x9,
        x10*x13*x19*x8-x10*x13*x20*x7-x10*x15*x17*x8+x10*x15*x20*x5+x10*x16*x17*x7-x10*x16*x19*x5-x11*x13*x18*x8+x11*x13*x20*x6+x11*x14*x17*x8-x11*x14*x20*x5-x11*x16*x17*x6+x11*x16*x18*x5+x12*x13*x18*x7-x12*x13*x19*x6-x12*x14*x17*x7+x12*x14*x19*x5+x12*x15*x17*x6-x12*x15*x18*x5-x14*x19*x8*x9+x14*x20*x7*x9+x15*x18*x8*x9-x15*x20*x6*x9-x16*x18*x7*x9+x16*x19*x6*x9,
        x1*x10*x15*x20-x1*x10*x16*x19-x1*x11*x14*x20+x1*x11*x16*x18+x1*x12*x14*x19-x1*x12*x15*x18+x10*x13*x19*x4-x10*x13*x20*x3-x10*x15*x17*x4+x10*x16*x17*x3-x11*x13*x18*x4+x11*x13*x2*x20+x11*x14*x17*x4-x11*x16*x17*x2+x12*x13*x18*x3-x12*x13*x19*x2-x12*x14*x17*x3+x12*x15*x17*x2-x14*x19*x4*x9+x14*x20*x3*x9+x15*x18*x4*x9-x15*x2*x20*x9-x16*x18*x3*x9+x16*x19*x2*x9,
        x1*x14*x19*x8-x1*x14*x20*x7-x1*x15*x18*x8+x1*x15*x20*x6+x1*x16*x18*x7-x1*x16*x19*x6+x13*x18*x3*x8-x13*x18*x4*x7-x13*x19*x2*x8+x13*x19*x4*x6+x13*x2*x20*x7-x13*x20*x3*x6-x14*x17*x3*x8+x14*x17*x4*x7-x14*x19*x4*x5+x14*x20*x3*x5+x15*x17*x2*x8-x15*x17*x4*x6+x15*x18*x4*x5-x15*x2*x20*x5-x16*x17*x2*x7+x16*x17*x3*x6-x16*x18*x3*x5+x16*x19*x2*x5,
        x1*x10*x19*x8-x1*x10*x20*x7-x1*x11*x18*x8+x1*x11*x20*x6+x1*x12*x18*x7-x1*x12*x19*x6-x10*x17*x3*x8+x10*x17*x4*x7-x10*x19*x4*x5+x10*x20*x3*x5+x11*x17*x2*x8-x11*x17*x4*x6+x11*x18*x4*x5-x11*x2*x20*x5-x12*x17*x2*x7+x12*x17*x3*x6-x12*x18*x3*x5+x12*x19*x2*x5+x18*x3*x8*x9-x18*x4*x7*x9-x19*x2*x8*x9+x19*x4*x6*x9+x2*x20*x7*x9-x20*x3*x6*x9,
        x1*x10*x15*x8-x1*x10*x16*x7-x1*x11*x14*x8+x1*x11*x16*x6+x1*x12*x14*x7-x1*x12*x15*x6-x10*x13*x3*x8+x10*x13*x4*x7-x10*x15*x4*x5+x10*x16*x3*x5+x11*x13*x2*x8-x11*x13*x4*x6+x11*x14*x4*x5-x11*x16*x2*x5-x12*x13*x2*x7+x12*x13*x3*x6-x12*x14*x3*x5+x12*x15*x2*x5+x14*x3*x8*x9-x14*x4*x7*x9-x15*x2*x8*x9+x15*x4*x6*x9+x16*x2*x7*x9-x16*x3*x6*x9
    ]
end

# Source: https://web.archive.org/web/20201202185136/http://www.cecm.sfu.ca/%7Erpearcea/mgb.html
#
# Critical pair degree sequence: 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18
function bayes148(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32) = np.polynomial_ring(
        k, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14, :x15, :x16, :x17, :x18, :x19, :x20, :x21, :x22, :x23, :x24, :x25, :x26, :x27, :x28, :x29, :x30, :x31, :x32], internal_ordering=internal_ordering)
    sys = [
        -x23*x32+x24*x31, -x22*x32+x24*x30, -x22*x31+x23*x30, -x21*x32+x24*x29, -x21*x31+x23*x29, -x21*x30+x22*x29, -x12*x32+x16*x28,  -x19*x28+x20*x27,
        -x11*x31+x15*x27, -x18*x28+x20*x26, -x18*x27+x19*x26, -x10*x30+x14*x26, -x17*x28+x20*x25, -x17*x27+x19*x25, -x17*x26+x18*x25, -x9*x29+x13*x25, x20*x8-x24*x4,
        -x17*x20-x17*x24-2*x17*x28-x17*x32+x18*x19+x18*x23+2*x18*x27+x18*x31+x19*x22+x19*x30-x20*x21-x20*x29-x21*x24-x21*x28-2*x21*x32+x22*x23+x22*x27+2*x22*x31+x23*x26-x24*x25-x25*x28-x25*x32+x26*x27+x26*x31+x27*x30-x28*x29-x29*x32+x30*x31,
        x19*x7-x23*x3, x18*x6-x2*x22, -x1*x21+x17*x5
    ]
end

# Source: https://web.archive.org/web/20201202185136/http://www.cecm.sfu.ca/%7Erpearcea/mgb.html
#
# Critical pair degree sequence: 7,6,7,8,9,8,8,8,9,8,8,8,8,8,8,8,8,9
function gametwo2(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    R, (p1,p2,p3,p4,p5,p6,p7) = np.polynomial_ring(k, [:p1,:p2,:p3,:p4,:p5,:p6,:p7], internal_ordering=internal_ordering)
    sys = [
        3821*p2*p3*p4*p5*p6*p7-7730*p2*p3*p4*p5*p6-164*p2*p3*p4*p5*p7-2536*p2*p3*p4*p6*p7-4321*p2*p3*p5*p6*p7-2161*p2*p4*p5*p6*p7-2188*p3*p4*p5*p6*p7-486*p2*p3*p4*p5+3491*p2*p3*p4*p6+4247*p2*p3*p5*p6+3528*p2*p4*p5*p6+2616*p3*p4*p5*p6-101*p2*p3*p4*p7+1765*p2*p3*p5*p7+258*p2*p4*p5*p7-378*p3*p4*p5*p7+1246*p2*p3*p6*p7+2320*p2*p4*p6*p7+1776*p3*p4*p6*p7+1715*p2*p5*p6*p7+728*p3*p5*p6*p7+842*p4*p5*p6*p7+69*p2*p3*p4-1660*p2*p3*p5+1863*p2*p4*p5+1520*p3*p4*p5-245*p2*p3*p6-804*p2*p4*p6-2552*p3*p4*p6-3152*p2*p5*p6+40*p3*p5*p6-1213*p4*p5*p6+270*p2*p3*p7-851*p2*p4*p7+327*p3*p4*p7-1151*p2*p5*p7+1035*p3*p5*p7-161*p4*p5*p7-230*p2*p6*p7-294*p3*p6*p7-973*p4*p6*p7-264*p5*p6*p7+874*p2*p3-2212*p2*p4+168*p3*p4+511*p2*p5-918*p3*p5-2017*p4*p5-76*p2*p6+465*p3*p6+1629*p4*p6+856*p5*p6-54*p2*p7-1355*p3*p7+227*p4*p7+77*p5*p7-220*p6*p7-696*p2+458*p3+486*p4+661*p5-650*p6+671*p7-439,
        -6157*p1*p3*p4*p5*p6*p7+13318*p1*p3*p4*p5*p6+5928*p1*p3*p4*p5*p7+1904*p1*p3*p4*p6*p7+2109*p1*p3*p5*p6*p7+8475*p1*p4*p5*p6*p7+2878*p3*p4*p5*p6*p7-8339*p1*p3*p4*p5-2800*p1*p3*p4*p6-9649*p1*p3*p5*p6-10964*p1*p4*p5*p6-4481*p3*p4*p5*p6+251*p1*p3*p4*p7-4245*p1*p3*p5*p7-7707*p1*p4*p5*p7-2448*p3*p4*p5*p7+1057*p1*p3*p6*p7-3605*p1*p4*p6*p7+546*p3*p4*p6*p7-3633*p1*p5*p6*p7-699*p3*p5*p6*p7-4126*p4*p5*p6*p7-730*p1*p3*p4+5519*p1*p3*p5+8168*p1*p4*p5+4366*p3*p4*p5+2847*p1*p3*p6+2058*p1*p4*p6-1416*p3*p4*p6+8004*p1*p5*p6+4740*p3*p5*p6+5361*p4*p5*p6-677*p1*p3*p7+1755*p1*p4*p7-760*p3*p4*p7+3384*p1*p5*p7+2038*p3*p5*p7+4119*p4*p5*p7+812*p1*p6*p7+11*p3*p6*p7+2022*p4*p6*p7+2642*p5*p6*p7+1276*p1*p3-1723*p1*p4+121*p3*p4-6456*p1*p5-3710*p3*p5-4525*p4*p5-2187*p1*p6-1559*p3*p6-848*p4*p6-4041*p5*p6-83*p1*p7-12*p3*p7-1180*p4*p7-2747*p5*p7-1970*p6*p7+2575*p1-161*p3+2149*p4+4294*p5+1687*p6+958*p7-1950,
        182*p1*p2*p4*p5*p6*p7-2824*p1*p2*p4*p5*p6-3513*p1*p2*p4*p5*p7-3386*p1*p2*p4*p6*p7-2330*p1*p2*p5*p6*p7-2838*p1*p4*p5*p6*p7+1294*p2*p4*p5*p6*p7+4764*p1*p2*p4*p5+1647*p1*p2*p4*p6+4221*p1*p2*p5*p6+814*p1*p4*p5*p6+2738*p2*p4*p5*p6+4057*p1*p2*p4*p7+2403*p1*p2*p5*p7+2552*p1*p4*p5*p7+471*p2*p4*p5*p7+448*p1*p2*p6*p7+2336*p1*p4*p6*p7+1617*p2*p4*p6*p7+2220*p1*p5*p6*p7-1543*p2*p5*p6*p7+402*p4*p5*p6*p7-5184*p1*p2*p4-3983*p1*p2*p5+44*p1*p4*p5-1327*p2*p4*p5-581*p1*p2*p6-389*p1*p4*p6-2722*p2*p4*p6+443*p1*p5*p6-2893*p2*p5*p6-154*p4*p5*p6-1277*p1*p2*p7-2018*p1*p4*p7-509*p2*p4*p7-1254*p1*p5*p7+602*p2*p5*p7-464*p4*p5*p7-647*p1*p6*p7+922*p2*p6*p7-1463*p4*p6*p7+729*p5*p6*p7+2665*p1*p2+591*p1*p4+981*p2*p4-444*p1*p5+1818*p2*p5-1985*p4*p5-1818*p1*p6+197*p2*p6+1038*p4*p6+340*p5*p6+399*p1*p7-835*p2*p7+787*p4*p7-753*p5*p7-221*p6*p7+481*p1+260*p2+1713*p4+1219*p5+794*p6+762*p7-1231,
        2923*p1*p2*p3*p5*p6*p7-4328*p1*p2*p3*p5*p6-3674*p1*p2*p3*p5*p7-3291*p1*p2*p3*p6*p7-4955*p1*p2*p5*p6*p7-8*p1*p3*p5*p6*p7-135*p2*p3*p5*p6*p7+7784*p1*p2*p3*p5+3471*p1*p2*p3*p6+1544*p1*p2*p5*p6+1607*p1*p3*p5*p6+1710*p2*p3*p5*p6+2434*p1*p2*p3*p7+1408*p1*p2*p5*p7-215*p1*p3*p5*p7+507*p2*p3*p5*p7+2208*p1*p2*p6*p7+1920*p1*p3*p6*p7-389*p2*p3*p6*p7+1304*p1*p5*p6*p7+2480*p2*p5*p6*p7+102*p3*p5*p6*p7-2683*p1*p2*p3-3508*p1*p2*p5-3505*p1*p3*p5-2400*p2*p3*p5-2236*p1*p2*p6-1727*p1*p3*p6-1354*p2*p3*p6-1022*p1*p5*p6-2099*p2*p5*p6-918*p3*p5*p6-495*p1*p2*p7-109*p1*p3*p7+474*p2*p3*p7+268*p1*p5*p7+1084*p2*p5*p7-190*p3*p5*p7-666*p1*p6*p7-497*p2*p6*p7+615*p3*p6*p7-912*p5*p6*p7+473*p1*p2+742*p1*p3+186*p2*p3+1021*p1*p5+2556*p2*p5+2312*p3*p5+1075*p1*p6+920*p2*p6+164*p3*p6+80*p5*p6-199*p1*p7-1270*p2*p7-1050*p3*p7-724*p5*p7+136*p6*p7+740*p1-474*p2+37*p3-1056*p5+303*p6+833*p7+736,
        4990*p1*p2*p3*p4*p6*p7-3067*p1*p2*p3*p4*p6-1661*p1*p2*p3*p4*p7-4064*p1*p2*p3*p6*p7-223*p1*p2*p4*p6*p7-5229*p1*p3*p4*p6*p7-4636*p2*p3*p4*p6*p7+5720*p1*p2*p3*p4+4872*p1*p2*p3*p6+1643*p1*p2*p4*p6+4536*p1*p3*p4*p6+2451*p2*p3*p4*p6+1264*p1*p2*p3*p7+70*p1*p2*p4*p7+2213*p1*p3*p4*p7+1734*p2*p3*p4*p7+1698*p1*p2*p6*p7+3799*p1*p3*p6*p7+1622*p2*p3*p6*p7+901*p1*p4*p6*p7-496*p2*p4*p6*p7+3782*p3*p4*p6*p7-5591*p1*p2*p3-1303*p1*p2*p4-6383*p1*p3*p4-2332*p2*p3*p4-3179*p1*p2*p6-6257*p1*p3*p6-3654*p2*p3*p6-1830*p1*p4*p6-1473*p2*p4*p6-3278*p3*p4*p6-1462*p1*p2*p7-1495*p1*p3*p7-468*p2*p3*p7-400*p1*p4*p7+431*p2*p4*p7-1907*p3*p4*p7-1547*p1*p6*p7-214*p2*p6*p7-1423*p3*p6*p7-1625*p4*p6*p7+5708*p1*p2+3809*p1*p3+2053*p2*p3+2824*p1*p4+1122*p2*p4+3653*p3*p4+3658*p1*p6+3001*p2*p6+3890*p3*p6+2371*p4*p6+602*p1*p7+185*p2*p7+899*p3*p7+963*p4*p7+560*p6*p7-4557*p1-3536*p2-1635*p3-2552*p4-2595*p6-207*p7+2740,
        -1407*p1*p2*p3*p4*p5*p7+4444*p1*p2*p3*p4*p5+2350*p1*p2*p3*p4*p7+5424*p1*p2*p3*p5*p7-2524*p1*p2*p4*p5*p7-985*p1*p3*p4*p5*p7-431*p2*p3*p4*p5*p7-2662*p1*p2*p3*p4-5342*p1*p2*p3*p5-39*p1*p2*p4*p5-2525*p1*p3*p4*p5-2650*p2*p3*p4*p5-3553*p1*p2*p3*p7-71*p1*p2*p4*p7-3268*p1*p3*p4*p7-1140*p2*p3*p4*p7-702*p1*p2*p5*p7-924*p1*p3*p5*p7-2198*p2*p3*p5*p7+4087*p1*p4*p5*p7+2709*p2*p4*p5*p7+587*p3*p4*p5*p7+968*p1*p2*p3-150*p1*p2*p4+909*p1*p3*p4+4587*p2*p3*p4+929*p1*p2*p5+1804*p1*p3*p5+2226*p2*p3*p5-916*p1*p4*p5+906*p2*p4*p5+2735*p3*p4*p5+1894*p1*p2*p7+2998*p1*p3*p7+1611*p2*p3*p7+304*p1*p4*p7-1601*p2*p4*p7+2066*p3*p4*p7-1971*p1*p5*p7-480*p2*p5*p7-500*p3*p5*p7-2617*p4*p5*p7-532*p1*p2+2016*p1*p3-2574*p2*p3+529*p1*p4-1251*p2*p4-2082*p3*p4+280*p1*p5-852*p2*p5-476*p3*p5-340*p4*p5-924*p1*p7+253*p2*p7-1090*p3*p7+170*p4*p7+1204*p5*p7-869*p1+1394*p2-264*p3+719*p4+219*p5-128*p7+506,
        -901*p1*p2*p3*p4*p5*p6+1805*p1*p2*p3*p4*p5-1103*p1*p2*p3*p4*p6-1746*p1*p2*p3*p5*p6-1968*p1*p2*p4*p5*p6+3957*p1*p3*p4*p5*p6+1293*p2*p3*p4*p5*p6-523*p1*p2*p3*p4-2498*p1*p2*p3*p5+693*p1*p2*p4*p5-2805*p1*p3*p4*p5-722*p2*p3*p4*p5-770*p1*p2*p3*p6+1088*p1*p2*p4*p6-232*p1*p3*p4*p6+2657*p2*p3*p4*p6+3281*p1*p2*p5*p6-1066*p1*p3*p5*p6+240*p2*p3*p5*p6-1174*p1*p4*p5*p6+1304*p2*p4*p5*p6-2070*p3*p4*p5*p6+2571*p1*p2*p3+115*p1*p2*p4+3899*p1*p3*p4-4641*p2*p3*p4-752*p1*p2*p5+1531*p1*p3*p5+1178*p2*p3*p5+11*p1*p4*p5-1144*p2*p4*p5-1701*p3*p4*p5+592*p1*p2*p6+1140*p1*p3*p6+130*p2*p3*p6+304*p1*p4*p6-2273*p2*p4*p6-1224*p3*p4*p6-2*p1*p5*p6-1090*p2*p5*p6+585*p3*p5*p6+670*p4*p5*p6-1867*p1*p2-4780*p1*p3+1079*p2*p3-2435*p1*p4+2901*p2*p4+2073*p3*p4+499*p1*p5+908*p2*p5+323*p3*p5+1631*p4*p5-966*p1*p6-315*p2*p6-481*p3*p6+759*p4*p6-595*p5*p6+3233*p1-1978*p2+729*p3-1184*p4-40*p5+446*p6+282
    ]
end

# Source: https://web.archive.org/web/20201202185136/http://www.cecm.sfu.ca/%7Erpearcea/mgb.html
function jason210(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    R, (x1, x2, x3, x4, x5, x6, x7, x8) = np.polynomial_ring(k, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8], internal_ordering=internal_ordering)
    sys = [
        x1^2*x3^4+x1*x2*x3^2*x5^2+x1*x2*x3*x4*x5*x7+x1*x2*x3*x4*x6*x8+x1*x2*x4^2*x6^2+x2^2*x4^4, x2^6, x1^6
    ]
end

# Source: https://web.archive.org/web/20201202185136/http://www.cecm.sfu.ca/%7Erpearcea/mgb.html
# Critical pair degree sequence: 3,4,5,6,7,7,8,8,7,7,7,8,8,8,8,9,10
function alea6(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    R, (x0, x1, x2, x3, x4, x5) = np.polynomial_ring(k, [:x0, :x1, :x2, :x3, :x4, :x5], internal_ordering=internal_ordering)
    sys = [
        5*x0^2*x3+37*x1*x3*x4+32*x1*x3*x5+21*x3*x5+55*x4*x5,
        39*x0*x1*x5+23*x1^2*x4+57*x1*x2*x4+56*x1*x4^2+10*x2^2+52*x3*x4*x5,
        33*x0^2*x3+51*x0^2+42*x0*x3*x5+51*x1^2*x4+32*x1*x3^2+x5^3,
        44*x0*x3^2+42*x1*x3+47*x1*x4^2+12*x2*x3+2*x2*x4*x5+43*x3*x4^2,
        49*x0^2*x2+11*x0*x1*x2+39*x0*x3*x4+44*x0*x3*x4+54*x0*x3+45*x1^2*x4,
        48*x0*x2*x3+2*x2^2*x3+59*x2^2*x5+17*x2+36*x3^3+45*x4
    ]
end

# Source: https://web.archive.org/web/20201202185136/http://www.cecm.sfu.ca/%7Erpearcea/mgb.html
#
# Critical pair degree sequence: 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33
function mayr42(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    R, (x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51) = np.polynomial_ring(
        k, [:x1,:x2,:x3,:x4,:x5,:x6,:x7,:x8,:x9,:x10,:x11,:x12,:x13,:x14,:x15,:x16,:x17,:x18,:x19,:x20,:x21,:x22,:x23,:x24,:x25,:x26,:x27,:x28,:x29,:x30,:x31,:x32,:x33,:x34,:x35,:x36,:x37,:x38,:x39,:x40,:x41,:x42,:x43,:x44,:x45,:x46,:x47,:x48,:x49,:x50,:x51], internal_ordering=internal_ordering)
    sys = [
        -x10*x51+x4*x49, x3*x48-x51*x9, x2*x47-x51*x8, x1*x46-x51*x7,
        x4*x44-x49*x9, x3*x43-x48*x8, x2*x42-x47*x7, x1*x41-x46*x6,
        x39*x4-x49*x9, x3*x38-x48*x8, x2*x37-x47*x7, x1*x36-x46*x6, x34*x9-x49*x9, x34*x4-x5*x51,
        x33*x8-x48*x8, x3*x33-x4*x51, x32*x7-x47*x7, x2*x32-x3*x51, x31*x6-x46*x6, x1*x31-x2*x51, 
        x14*x39*x9-x29*x44*x9, x13*x38*x8-x28*x43*x8, x12*x37*x7-x27*x42*x7, x11*x36*x6-x26*x41*x6,
        x26^2*x46*x6-x51^3*x7, x11^2*x46*x6-x2*x51^3, x21^2*x41*x6-x46*x51^2*x6, x16^2*x36*x6-x46*x51^2*x6,
        x24*x30*x39*x50*x9-x29*x44*x50*x51*x9, x23*x29*x38*x49*x8-x28*x43*x49*x51*x8, x22*x28*x37*x48*x7-x27*x42*x48*x51*x7,
        x21*x27*x36*x47*x6-x26*x41*x47*x51*x6, x24*x25*x39*x45*x9-x29*x44*x45*x51*x9, x23*x24*x38*x44*x8-x28*x43*x44*x51*x8,
        x22*x23*x37*x43*x7-x27*x42*x43*x51*x7, x21*x22*x36*x42*x6-x26*x41*x42*x51*x6, x20*x24*x39*x40*x9-x29*x40*x44*x51*x9,
        x19*x23*x38*x39*x8-x28*x39*x43*x51*x8, x15*x24*x35*x39*x9-x29*x35*x44*x51*x9, x18*x22*x37*x38*x7-x27*x38*x42*x51*x7, 
        x14*x23*x34*x38*x8-x28*x34*x43*x51*x8, x17*x21*x36*x37*x6-x26*x37*x41*x51*x6, x13*x22*x33*x37*x7-x27*x33*x42*x51*x7,
        x12*x21*x32*x36*x6-x26*x32*x41*x51*x6
    ]
end

###
# Some other examples

function ku10(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) =
        np.polynomial_ring(k, ["x$i" for i in 1:10], internal_ordering=internal_ordering)
    [
        5 * x1 * x2 + 5 * x1 + 3 * x2 + 55,
        7 * x2 * x3 + 9 * x2 + 9 * x3 + 19,
        3 * x3 * x4 + 6 * x3 + 5 * x4 - 4,
        6 * x4 * x5 + 6 * x4 + 7 * x5 + 118,
        x5 * x6 + 3 * x5 + 9 * x6 + 27,
        6 * x6 * x7 + 7 * x6 + x7 + 72,
        9 * x7 * x8 + 7 * x7 + x8 + 35,
        4 * x8 * x9 + 4 * x8 + 6 * x9 + 16,
        8 * x9 * x10 + 4 * x9 + 3 * x10 - 51,
        3 * x1 * x10 - 6 * x1 + x10 + 5
    ]
end

function kinema(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (z1, z2, z3, z4, z5, z6, z7, z8, z9) =
        np.polynomial_ring(k, ["z$i" for i in 1:9], internal_ordering=internal_ordering)
    [
        z1^2 + z2^2 + z3^2 - 12 * z1 - 68
        z4^2 + z5^2 + z6^2 - 12 * z5 - 68
        z7^2 + z8^2 + z9^2 - 24 * z8 - 12 * z9 + 100
        z1 * z4 + z2 * z5 + z3 * z6 - 6 * z1 - 6 * z5 - 52
        z1 * z7 + z2 * z8 + z3 * z9 - 6 * z1 - 12 * z8 - 6 * z9 + 64
        z4 * z7 + z5 * z8 + z6 * z9 - 6 * z5 - 12 * z8 - 6 * z9 + 32
        2 * z2 + 2 * z3 - z4 - z5 - 2 * z6 - z7 - z9 + 18
        z1 + z2 + 2 * z3 + 2 * z4 + 2 * z6 - 2 * z7 + z8 - z9 - 38
        z1 + z3 - 2 * z4 + z5 - z6 + 2 * z7 - 2 * z8 + 8
    ]
end

function s9_1(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (a, b, c, d, e, f, g, h) =
        np.polynomial_ring(k, ["x$i" for i in 1:8], internal_ordering=internal_ordering)
    [
        -e * g - 2 * d * h,
        9 * e + 4 * b,
        -4 * c * h - 2 * e * f - 3 * d * g,
        -7 * c + 9 * a - 8 * f,
        -4 * d * f - 5 * c * g - 6 * h - 3 * e,
        -5 * d - 6 * c * f - 7 * g + 9 * b,
        9 * d + 6 * a - 5 * b,
        9 * c - 7 * a + 8
    ]
end

function ojika4(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x1, x2, x3) = np.polynomial_ring(k, ["x$i" for i in 1:3], internal_ordering=internal_ordering)
    [
        x1 + x3 * x1^3 + x1 * x3 * x2^2 - x1 * x3,
        10 * x2 - 2 * x2 * x3 * x1^2 - x3 * x2^3 - x2 * x3,
        -6 * x3^2 * x1^4 - 3 * x1^2 * x2^2 * x3^2 - x3^2 * x1^2 + 28 * x3 * x1^2 -
        3 * x3^2 * x2^4 +
        2 * x3^2 * x2^2 +
        7 * x3 * x2^2 +
        x3^2 - 11 * x3 + 10
    ]
end

function ojika3_d1R2(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x1, x2, x3) = np.polynomial_ring(k, ["x$i" for i in 1:3], internal_ordering=internal_ordering)
    [
        x1^3 * x3 + x1 * x3 * x2^2 - x1 * x3 + x1,
        -2 * x1^2 * x3 * x2 - x3 * x2^3 - x3 * x2 + 10 * x2,
        -6 * x1^4 * x3^2 - 3 * x1^2 * x3^2 * x2^2 - 3 * x3^2 * x2^4 - x1^2 * x3^2 +
        2 * x3^2 * x2^2 +
        28 * x1^2 * x3 +
        7 * x3 * x2^2 +
        x3^2 - 11 * x3 + 10
    ]
end

function ojika4_d1R2_d2R5(; np=AbstractAlgebra, k=np.QQ, internal_ordering=:degrevlex)
    _, (x1, x2, x3) = np.polynomial_ring(k, ["x$i" for i in 1:3], internal_ordering=internal_ordering)
    [
        x1^3 * x3 + x1 * x3 * x2^2 - x1 * x3 + x1,
        -2 * x1^2 * x3 * x2 - x3 * x2^3 - x3 * x2 + 10 * x2,
        -6 * x1^4 * x3^2 - 3 * x1^2 * x3^2 * x2^2 - 3 * x3^2 * x2^4 - x1^2 * x3^2 +
        2 * x3^2 * x2^2 +
        28 * x1^2 * x3 +
        7 * x3 * x2^2 +
        x3^2 - 11 * x3 + 10
    ]
end

# SIAN examples

# Critical pair degree sequence: 3,3,4,3,2,2,2,3,3,2,4,4,5,4,3,4,3,4,4,4,3,3,2,3,2,2,3,3,4,3,2,3,5,6,7
function HIV2(; np=AbstractAlgebra, internal_ordering=:degrevlex, k=np.QQ)
    R, (z_9,w_8,y_8,z_8,x_7,v_7,w_7,y_7,z_7,x_6,v_6,w_6,y_6,z_6,x_5,v_5,w_5,y_5,z_5,x_4,v_4,w_4,y_4,z_4,x_3,v_3,w_3,y_3,z_3,x_2,v_2,w_2,y_2,z_2,x_1,v_1,w_1,y_1,z_1,x_0,v_0,w_0,y_0,z_0,z_aux,lm_0,d_0,beta_0,a_0,k_0,u_0,c_0,q_0,b_0,h_0) = np.polynomial_ring(k, [:z_9,:w_8,:y_8,:z_8,:x_7,:v_7,:w_7,:y_7,:z_7,:x_6,:v_6,:w_6,:y_6,:z_6,:x_5,:v_5,:w_5,:y_5,:z_5,:x_4,:v_4,:w_4,:y_4,:z_4,:x_3,:v_3,:w_3,:y_3,:z_3,:x_2,:v_2,:w_2,:y_2,:z_2,:x_1,:v_1,:w_1,:y_1,:z_1,:x_0,:v_0,:w_0,:y_0,:z_0,:z_aux,:lm_0,:d_0,:beta_0,:a_0,:k_0,:u_0,:c_0,:q_0,:b_0,:h_0], internal_ordering=internal_ordering)
    sys = [
    		-z_0 + 391579310130776933615871,
		-w_0*y_0*c_0*q_0 + z_0*h_0 + z_1,
		-w_0 + 11747232036185701625123157,
		-w_0*y_0*z_0*c_0 + w_0*y_0*c_0*q_0 + w_0*b_0 + w_1,
		-z_1 + 603238345346139528359522036104641908547430757404926535076096396330351543065342569329178192530970096,
		-y_1*w_0*c_0*q_0 - w_1*y_0*c_0*q_0 + z_1*h_0 + z_2,
		-x_0*v_0*beta_0 + y_0*a_0 + y_1,
		-w_1 - 574390083393820420103212759171522849275323800567457691692315556508774738634619943653674738585986743,
		-z_1*w_0*y_0*c_0 - y_1*w_0*z_0*c_0 - w_1*y_0*z_0*c_0 + y_1*w_0*c_0*q_0 + w_1*y_0*c_0*q_0 + w_1*b_0 + w_2,
		-z_2 - 29495809942494772235717487304183887142454587420869561152232261111077119129926946003670610574309564618092345101563794158451252951899401078528644667279465056622989870416543178,
		-2*w_1*y_1*c_0*q_0 - y_2*w_0*c_0*q_0 - w_2*y_0*c_0*q_0 + z_2*h_0 + z_3,
		-v_1*x_0*beta_0 - x_1*v_0*beta_0 + y_1*a_0 + y_2,
		-y_0*k_0 + v_0*u_0 + v_1,
		x_0*v_0*beta_0 + x_0*d_0 + x_1 - lm_0,
		-w_2 + 72526766702336286507825226434003582721065918205231512442169043118987526401873302266876059227470399469758729806573577534508525039097990631413723803846283476257933579221333353,
		-2*y_1*z_1*w_0*c_0 - 2*w_1*z_1*y_0*c_0 - z_2*w_0*y_0*c_0 - 2*w_1*y_1*z_0*c_0 - y_2*w_0*z_0*c_0 - w_2*y_0*z_0*c_0 + 2*w_1*y_1*c_0*q_0 + y_2*w_0*c_0*q_0 + w_2*y_0*c_0*q_0 + w_2*b_0 + w_3,
		-z_3 + 3724360479477568236622724779173882605002730995788446608513197484643942968717402455098831776126940201242074854832830505138391424197393570926108499423090983113588950610968169165079829434421433538510808865018487980237424646729309138096329558471260151,
		-3*y_2*w_1*c_0*q_0 - 3*w_2*y_1*c_0*q_0 - y_3*w_0*c_0*q_0 - w_3*y_0*c_0*q_0 + z_3*h_0 + z_4,
		-2*x_1*v_1*beta_0 - v_2*x_0*beta_0 - x_2*v_0*beta_0 + y_2*a_0 + y_3,
		-y_1*k_0 + v_1*u_0 + v_2,
		v_1*x_0*beta_0 + x_1*v_0*beta_0 + x_1*d_0 + x_2,
		-w_3 - 10065260635460367136827570054743039971303285538244550872182162715368975218110847370878710258279443731462271690633257400089986872117350319845914336108436246847573437942215246570949277640001082095423119390494264184527042153943726452445744014647482237,
		-6*w_1*y_1*z_1*c_0 - 3*z_2*y_1*w_0*c_0 - 3*y_2*z_1*w_0*c_0 - 3*z_2*w_1*y_0*c_0 - 3*w_2*z_1*y_0*c_0 - z_3*w_0*y_0*c_0 - 3*y_2*w_1*z_0*c_0 - 3*w_2*y_1*z_0*c_0 - y_3*w_0*z_0*c_0 - w_3*y_0*z_0*c_0 + 3*y_2*w_1*c_0*q_0 + 3*w_2*y_1*c_0*q_0 + y_3*w_0*c_0*q_0 + w_3*y_0*c_0*q_0 + w_3*b_0 + w_4,
		-z_4 - 516866539497097390741430822517376417656042385836133479039598482854041835920781599692935765220288317208568307674240570033041820302571872340249173134476132079789328276084646776952006251103233801864353824418884420269152114838434532519317164510254202954894504821067054540160787808587306789965526678376390381380646541732318163,
		-6*w_2*y_2*c_0*q_0 - 4*y_3*w_1*c_0*q_0 - 4*w_3*y_1*c_0*q_0 - y_4*w_0*c_0*q_0 - w_4*y_0*c_0*q_0 + z_4*h_0 + z_5,
		-3*v_2*x_1*beta_0 - 3*x_2*v_1*beta_0 - v_3*x_0*beta_0 - x_3*v_0*beta_0 + y_3*a_0 + y_4,
		2*x_1*v_1*beta_0 + v_2*x_0*beta_0 + x_2*v_0*beta_0 + x_2*d_0 + x_3,
		-y_2*k_0 + v_2*u_0 + v_3,
		-z_5 + 98000219445755281465191234047409519457731340099562049701122536693667933661439374149287916104376587314997709305518514997881274053044896236528087300963711221079108242783818349292217649914304832608846704663912984082981256664618016281150945974695812571637338872669364928091333910816303015480587200779217401127972575362198750665563610538189010755829185777812616412129948629004723126787453477990542388,
		-10*y_3*w_2*c_0*q_0 - 10*w_3*y_2*c_0*q_0 - 5*y_4*w_1*c_0*q_0 - 5*w_4*y_1*c_0*q_0 - y_5*w_0*c_0*q_0 - w_5*y_0*c_0*q_0 + z_5*h_0 + z_6,
		-6*x_2*v_2*beta_0 - 4*v_3*x_1*beta_0 - 4*x_3*v_1*beta_0 - v_4*x_0*beta_0 - x_4*v_0*beta_0 + y_4*a_0 + y_5,
		-12*z_2*w_1*y_1*c_0 - 12*y_2*w_1*z_1*c_0 - 12*w_2*y_1*z_1*c_0 - 6*y_2*z_2*w_0*c_0 - 4*z_3*y_1*w_0*c_0 - 4*y_3*z_1*w_0*c_0 - 6*w_2*z_2*y_0*c_0 - 4*z_3*w_1*y_0*c_0 - 4*w_3*z_1*y_0*c_0 - z_4*w_0*y_0*c_0 - 6*w_2*y_2*z_0*c_0 - 4*y_3*w_1*z_0*c_0 - 4*w_3*y_1*z_0*c_0 - y_4*w_0*z_0*c_0 - w_4*y_0*z_0*c_0 + 6*w_2*y_2*c_0*q_0 + 4*y_3*w_1*c_0*q_0 + 4*w_3*y_1*c_0*q_0 + y_4*w_0*c_0*q_0 + w_4*y_0*c_0*q_0 + w_4*b_0 + w_5,
		3*v_2*x_1*beta_0 + 3*x_2*v_1*beta_0 + v_3*x_0*beta_0 + x_3*v_0*beta_0 + x_3*d_0 + x_4,
		-y_3*k_0 + v_3*u_0 + v_4,
		-z_6 - 21458030562658662580738258553823341766587663576114394834695028559774918019820596509904405389041607988467364762358827103883264176228057170978494313118084991093862092910419070546968589068608627622285541608722719273109079816490152719189187972967958945957946933218486487309661517029087886222602047679212036481982516143941270331116995678709894046484339971261780233952395479938807125800855113714268132010322034207128699220553082030907388119389480410220733479809310348610218079,
		-20*w_3*y_3*c_0*q_0 - 15*y_4*w_2*c_0*q_0 - 15*w_4*y_2*c_0*q_0 - 6*y_5*w_1*c_0*q_0 - 6*w_5*y_1*c_0*q_0 - y_6*w_0*c_0*q_0 - w_6*y_0*c_0*q_0 + z_6*h_0 + z_7,
		-10*v_3*x_2*beta_0 - 10*x_3*v_2*beta_0 - 5*v_4*x_1*beta_0 - 5*x_4*v_1*beta_0 - v_5*x_0*beta_0 - x_5*v_0*beta_0 + y_5*a_0 + y_6,
		-30*y_2*z_2*w_1*c_0 - 30*w_2*z_2*y_1*c_0 - 20*z_3*w_1*y_1*c_0 - 30*w_2*y_2*z_1*c_0 - 20*y_3*w_1*z_1*c_0 - 20*w_3*y_1*z_1*c_0 - 10*z_3*y_2*w_0*c_0 - 10*y_3*z_2*w_0*c_0 - 5*z_4*y_1*w_0*c_0 - 5*y_4*z_1*w_0*c_0 - 10*z_3*w_2*y_0*c_0 - 10*w_3*z_2*y_0*c_0 - 5*z_4*w_1*y_0*c_0 - 5*w_4*z_1*y_0*c_0 - z_5*w_0*y_0*c_0 - 10*y_3*w_2*z_0*c_0 - 10*w_3*y_2*z_0*c_0 - 5*y_4*w_1*z_0*c_0 - 5*w_4*y_1*z_0*c_0 - y_5*w_0*z_0*c_0 - w_5*y_0*z_0*c_0 + 10*y_3*w_2*c_0*q_0 + 10*w_3*y_2*c_0*q_0 + 5*y_4*w_1*c_0*q_0 + 5*w_4*y_1*c_0*q_0 + y_5*w_0*c_0*q_0 + w_5*y_0*c_0*q_0 + w_5*b_0 + w_6,
		6*x_2*v_2*beta_0 + 4*v_3*x_1*beta_0 + 4*x_3*v_1*beta_0 + v_4*x_0*beta_0 + x_4*v_0*beta_0 + x_4*d_0 + x_5,
		-y_4*k_0 + v_4*u_0 + v_5,
		-z_7 + 5577749280924487080075757231171031743354822071715892112217952137532438336218771414445727554974797200816590802671571482512393969909942909236616088107078369196678728899070868401051382034652738349436714476525941875879020437352974716904364501936790077254158995922237738283579804264268588965519025664074896789072276888776123306216048524032102523476927992446263909223919108230846122070777164234943059560088723172330230859379939527427325462452530168613108231688004536190574621583011662698450824506643742615018133550422437856474198827180129780321053982,
		-35*y_4*w_3*c_0*q_0 - 35*w_4*y_3*c_0*q_0 - 21*y_5*w_2*c_0*q_0 - 21*w_5*y_2*c_0*q_0 - 7*y_6*w_1*c_0*q_0 - 7*w_6*y_1*c_0*q_0 - y_7*w_0*c_0*q_0 - w_7*y_0*c_0*q_0 + z_7*h_0 + z_8,
		-90*w_2*y_2*z_2*c_0 - 60*z_3*y_2*w_1*c_0 - 60*y_3*z_2*w_1*c_0 - 60*z_3*w_2*y_1*c_0 - 60*w_3*z_2*y_1*c_0 - 30*z_4*w_1*y_1*c_0 - 60*y_3*w_2*z_1*c_0 - 60*w_3*y_2*z_1*c_0 - 30*y_4*w_1*z_1*c_0 - 30*w_4*y_1*z_1*c_0 - 20*y_3*z_3*w_0*c_0 - 15*z_4*y_2*w_0*c_0 - 15*y_4*z_2*w_0*c_0 - 6*z_5*y_1*w_0*c_0 - 6*y_5*z_1*w_0*c_0 - 20*w_3*z_3*y_0*c_0 - 15*z_4*w_2*y_0*c_0 - 15*w_4*z_2*y_0*c_0 - 6*z_5*w_1*y_0*c_0 - 6*w_5*z_1*y_0*c_0 - z_6*w_0*y_0*c_0 - 20*w_3*y_3*z_0*c_0 - 15*y_4*w_2*z_0*c_0 - 15*w_4*y_2*z_0*c_0 - 6*y_5*w_1*z_0*c_0 - 6*w_5*y_1*z_0*c_0 - y_6*w_0*z_0*c_0 - w_6*y_0*z_0*c_0 + 20*w_3*y_3*c_0*q_0 + 15*y_4*w_2*c_0*q_0 + 15*w_4*y_2*c_0*q_0 + 6*y_5*w_1*c_0*q_0 + 6*w_5*y_1*c_0*q_0 + y_6*w_0*c_0*q_0 + w_6*y_0*c_0*q_0 + w_6*b_0 + w_7,
		-20*x_3*v_3*beta_0 - 15*v_4*x_2*beta_0 - 15*x_4*v_2*beta_0 - 6*v_5*x_1*beta_0 - 6*x_5*v_1*beta_0 - v_6*x_0*beta_0 - x_6*v_0*beta_0 + y_6*a_0 + y_7,
		10*v_3*x_2*beta_0 + 10*x_3*v_2*beta_0 + 5*v_4*x_1*beta_0 + 5*x_4*v_1*beta_0 + v_5*x_0*beta_0 + x_5*v_0*beta_0 + x_5*d_0 + x_6,
		-y_5*k_0 + v_5*u_0 + v_6,
		-z_8 - 1644205192807033404764559007893172745748849984384312802354200787599776059355428090998225620786645231842428680599803099279698065002526758747161431644262518692015208453655878780813913561053608242320460813887802099921967392102887109597652291416944346692768469655604042964417877143328586109034623413791328007722399916097968656202908830017497359716003032873011341768633958462210430990403673940581962497822841126154764133423680404948778245115721750781412392929302931078597797754293988416603816534949548588565217177356266284435209162765551819170775850481304775023438958310814160758575680783396565217760341139831156750260898232,
		-70*w_4*y_4*c_0*q_0 - 56*y_5*w_3*c_0*q_0 - 56*w_5*y_3*c_0*q_0 - 28*y_6*w_2*c_0*q_0 - 28*w_6*y_2*c_0*q_0 - 8*y_7*w_1*c_0*q_0 - 8*w_7*y_1*c_0*q_0 - y_8*w_0*c_0*q_0 - w_8*y_0*c_0*q_0 + z_8*h_0 + z_9,
		-35*v_4*x_3*beta_0 - 35*x_4*v_3*beta_0 - 21*v_5*x_2*beta_0 - 21*x_5*v_2*beta_0 - 7*v_6*x_1*beta_0 - 7*x_6*v_1*beta_0 - v_7*x_0*beta_0 - x_7*v_0*beta_0 + y_7*a_0 + y_8,
		-210*z_3*w_2*y_2*c_0 - 210*y_3*w_2*z_2*c_0 - 210*w_3*y_2*z_2*c_0 - 140*y_3*z_3*w_1*c_0 - 105*z_4*y_2*w_1*c_0 - 105*y_4*z_2*w_1*c_0 - 140*w_3*z_3*y_1*c_0 - 105*z_4*w_2*y_1*c_0 - 105*w_4*z_2*y_1*c_0 - 42*z_5*w_1*y_1*c_0 - 140*w_3*y_3*z_1*c_0 - 105*y_4*w_2*z_1*c_0 - 105*w_4*y_2*z_1*c_0 - 42*y_5*w_1*z_1*c_0 - 42*w_5*y_1*z_1*c_0 - 35*z_4*y_3*w_0*c_0 - 35*y_4*z_3*w_0*c_0 - 21*z_5*y_2*w_0*c_0 - 21*y_5*z_2*w_0*c_0 - 7*z_6*y_1*w_0*c_0 - 7*y_6*z_1*w_0*c_0 - 35*z_4*w_3*y_0*c_0 - 35*w_4*z_3*y_0*c_0 - 21*z_5*w_2*y_0*c_0 - 21*w_5*z_2*y_0*c_0 - 7*z_6*w_1*y_0*c_0 - 7*w_6*z_1*y_0*c_0 - z_7*w_0*y_0*c_0 - 35*y_4*w_3*z_0*c_0 - 35*w_4*y_3*z_0*c_0 - 21*y_5*w_2*z_0*c_0 - 21*w_5*y_2*z_0*c_0 - 7*y_6*w_1*z_0*c_0 - 7*w_6*y_1*z_0*c_0 - y_7*w_0*z_0*c_0 - w_7*y_0*z_0*c_0 + 35*y_4*w_3*c_0*q_0 + 35*w_4*y_3*c_0*q_0 + 21*y_5*w_2*c_0*q_0 + 21*w_5*y_2*c_0*q_0 + 7*y_6*w_1*c_0*q_0 + 7*w_6*y_1*c_0*q_0 + y_7*w_0*c_0*q_0 + w_7*y_0*c_0*q_0 + w_7*b_0 + w_8,
		-y_6*k_0 + v_6*u_0 + v_7,
		20*x_3*v_3*beta_0 + 15*v_4*x_2*beta_0 + 15*x_4*v_2*beta_0 + 6*v_5*x_1*beta_0 + 6*x_5*v_1*beta_0 + v_6*x_0*beta_0 + x_6*v_0*beta_0 + x_6*d_0 + x_7,
		-z_9 + 547169434379018145028611952347265787172199439668093986338563317785852393916605948818072737647734930569662248686887630106250170916620431557805104040465245886891859630463989194155254027609865997719559075910769153141887632327231757960445334365558978268855512815198121035198224568168185312239960245138751315716963041810932606675513939807949557458205653521468787914237949369491510319720103583796426594093631148917030950116046047765636929109396208199320896283070375401931791154858990255781028687359602238037994809564110293220500908246017866610361218898309281015955618403398237167794184297295140520001201453021971553871890098816334741387420433906924143868197736158620363296005190522485295124049435107,
		-w_4 + 1908418664542663099774091716482431588374718871839089056219614579971723103702336041333144915382378846941799838921829645861013275662853484016597679102928009974796834453662803702825468863863699315044327414900631216745549361241879675632162230381413209813079885090508303150179123212793333319289763289097064792543950982117615554,
		-w_5 - 417865452360254033504908453978493737666079490808113004239654956627219414784331711209111925283255385975513322275446603895846428020148233461120142782852541875605278407949821612766657667057555188148906991687680366739791207569112150140359955591390126222672894978054822180835696002982092410893734533201028323111216374274942030006580232521518488130704516322548470273619123822892391160346504020692095957,
		-w_6 + 108618948958045071996686191814880701723273779780900537851211480014575954859862235158064614208312594040667185441627007847229742411532034623787664180692273256112833274203152565248107495837749917833198273389140944246934942898473699518366663922813793434057738114266639202852428237429743321235897912960657393304612016304368015117310995538286830206980622039291214698553695303384645857319518546092967736770353565298943410835213622034668330539695564620041933350395927208045146375,
		-w_7 - 32018620938167915709424800884179286400769850246168416971678757387326577790632020350923544485495646077637968187175721185871067568698114797069132849683581823639272546363658830702006575242653103871085066129287516159880704931949640223431668861822555144067653403459119184691919550430168307856142339035519933356454825520213063466807127469721496469647277108997670176435798508951676198010280591284129007081055926705383684593216547231892829664223470094372139619480976098074519244480072648434009312814289977866073070383674337616936127106839527611963606981,
		-w_8 + 10655367581234525241629788927939090277671296368705864939130569241961936206395813457389337876591554076728476295485298661373225491120223759292480329085090132276945483193842077671935435551158231403512950548263059434981109528221157431942640281970058313649122539084766043649565177284923008382695404456058860264634269420285791899882734001989063871991771314320957623245853109055586163026501019523045737918594445219682828379192621279639096563176615776544296906282357060466021309911440499438196138881527087217758263306505943656173773494467108160684251721889266292373387780029621670660070518810549666699067877947310419393249067624,
		z_aux - 1
    ]
end

# Critical pair degree sequence: 3,2,3,3,3,3,3,3,3,3,4,2,3,4,5
function ChemicalReactionNetwork(; np=AbstractAlgebra, internal_ordering=:degrevlex, k=np.QQ)
    R, (x2_6,x3_6,x4_5,x2_5,x5_5,x1_5,x3_5,x6_5,x4_4,x2_4,x5_4,x1_4,x3_4,x6_4,x4_3,x2_3,x5_3,x1_3,x3_3,x6_3,x4_2,x2_2,x5_2,x1_2,x3_2,x6_2,x4_1,x2_1,x5_1,x1_1,x3_1,x6_1,x4_0,x2_0,x5_0,x1_0,x3_0,x6_0,z_aux,k1_0,k2_0,k4_0,k3_0,k5_0,k6_0) = np.polynomial_ring(k, [:x2_6,:x3_6,:x4_5,:x2_5,:x5_5,:x1_5,:x3_5,:x6_5,:x4_4,:x2_4,:x5_4,:x1_4,:x3_4,:x6_4,:x4_3,:x2_3,:x5_3,:x1_3,:x3_3,:x6_3,:x4_2,:x2_2,:x5_2,:x1_2,:x3_2,:x6_2,:x4_1,:x2_1,:x5_1,:x1_1,:x3_1,:x6_1,:x4_0,:x2_0,:x5_0,:x1_0,:x3_0,:x6_0,:z_aux,:k1_0,:k2_0,:k4_0,:k3_0,:k5_0,:k6_0], internal_ordering=internal_ordering)
    sys = [
    		-x2_0 + 544397587489487158925,
		-x2_0*x1_0*k1_0 - x4_0*k2_0 - x4_0*k3_0 + x2_1,
		-x3_0 + 29682033801657517034,
		x5_0*x3_0*k6_0 - x4_0*k3_0 - x6_0*k5_0 + x3_1,
		-x2_1 + 31319937365291324106855334029149083433729222633545127704567987,
		-x1_1*x2_0*k1_0 - x2_1*x1_0*k1_0 - x4_1*k2_0 - x4_1*k3_0 + x2_2,
		-x2_0*x1_0*k1_0 + x4_0*k2_0 + x4_0*k3_0 + x4_1,
		x2_0*x1_0*k1_0 - x4_0*k2_0 - x6_0*k4_0 + x1_1,
		-x3_1 - 1205403898926015103232327521221029715967312994127374687044129,
		x3_1*x5_0*k6_0 + x5_1*x3_0*k6_0 - x4_1*k3_0 - x6_1*k5_0 + x3_2,
		-x5_0*x3_0*k6_0 + x6_0*k4_0 + x6_0*k5_0 + x6_1,
		x5_0*x3_0*k6_0 - x6_0*k4_0 - x6_0*k5_0 + x5_1,
		-x2_2 - 229339524796362821697179199520342673234123591683553488185657952768354303216823322458483803211712707633,
		-2*x2_1*x1_1*k1_0 - x1_2*x2_0*k1_0 - x2_2*x1_0*k1_0 - x4_2*k2_0 - x4_2*k3_0 + x2_3,
		-x1_1*x2_0*k1_0 - x2_1*x1_0*k1_0 + x4_1*k2_0 + x4_1*k3_0 + x4_2,
		x1_1*x2_0*k1_0 + x2_1*x1_0*k1_0 - x4_1*k2_0 - x6_1*k4_0 + x1_2,
		-x3_2 + 51994373568763482184084092961808697028251976190716290042832706767730767453744791804300393099736005865,
		2*x5_1*x3_1*k6_0 + x3_2*x5_0*k6_0 + x5_2*x3_0*k6_0 - x4_2*k3_0 - x6_2*k5_0 + x3_3,
		x3_1*x5_0*k6_0 + x5_1*x3_0*k6_0 - x6_1*k4_0 - x6_1*k5_0 + x5_2,
		-x3_1*x5_0*k6_0 - x5_1*x3_0*k6_0 + x6_1*k4_0 + x6_1*k5_0 + x6_2,
		-x2_3 - 232038201156430746537091630904147281723828787963263084910508601531237857942582213788880858182009072907809417201406568624315583493857533204297241,
		-3*x1_2*x2_1*k1_0 - 3*x2_2*x1_1*k1_0 - x1_3*x2_0*k1_0 - x2_3*x1_0*k1_0 - x4_3*k2_0 - x4_3*k3_0 + x2_4,
		2*x2_1*x1_1*k1_0 + x1_2*x2_0*k1_0 + x2_2*x1_0*k1_0 - x4_2*k2_0 - x6_2*k4_0 + x1_3,
		-2*x2_1*x1_1*k1_0 - x1_2*x2_0*k1_0 - x2_2*x1_0*k1_0 + x4_2*k2_0 + x4_2*k3_0 + x4_3,
		-x3_3 - 2489841093515072908413476041338707895889602410462016515249818481973712038712671501796124969357657634299055404533594664232527151976813331542405,
		3*x3_2*x5_1*k6_0 + 3*x5_2*x3_1*k6_0 + x3_3*x5_0*k6_0 + x5_3*x3_0*k6_0 - x4_3*k3_0 - x6_3*k5_0 + x3_4,
		2*x5_1*x3_1*k6_0 + x3_2*x5_0*k6_0 + x5_2*x3_0*k6_0 - x6_2*k4_0 - x6_2*k5_0 + x5_3,
		-2*x5_1*x3_1*k6_0 - x3_2*x5_0*k6_0 - x5_2*x3_0*k6_0 + x6_2*k4_0 + x6_2*k5_0 + x6_3,
		-x2_4 + 6833268318610526294398509674052602741963506552428132533026349353865937016318327644539867625341607863425291998178353893918539978698743728647984446868290545116877979367133284448342699619,
		-6*x2_2*x1_2*k1_0 - 4*x1_3*x2_1*k1_0 - 4*x2_3*x1_1*k1_0 - x1_4*x2_0*k1_0 - x2_4*x1_0*k1_0 - x4_4*k2_0 - x4_4*k3_0 + x2_5,
		3*x1_2*x2_1*k1_0 + 3*x2_2*x1_1*k1_0 + x1_3*x2_0*k1_0 + x2_3*x1_0*k1_0 - x4_3*k2_0 - x6_3*k4_0 + x1_4,
		-3*x1_2*x2_1*k1_0 - 3*x2_2*x1_1*k1_0 - x1_3*x2_0*k1_0 - x2_3*x1_0*k1_0 + x4_3*k2_0 + x4_3*k3_0 + x4_4,
		-x3_4 + 139372684659646947400032617915256822925543857095425418659231178320203702886011430928972731874199419309231800473513836491566273073235642297304567870916637554479092119674738526118464349,
		6*x5_2*x3_2*k6_0 + 4*x3_3*x5_1*k6_0 + 4*x5_3*x3_1*k6_0 + x3_4*x5_0*k6_0 + x5_4*x3_0*k6_0 - x4_4*k3_0 - x6_4*k5_0 + x3_5,
		3*x3_2*x5_1*k6_0 + 3*x5_2*x3_1*k6_0 + x3_3*x5_0*k6_0 + x5_3*x3_0*k6_0 - x6_3*k4_0 - x6_3*k5_0 + x5_4,
		-3*x3_2*x5_1*k6_0 - 3*x5_2*x3_1*k6_0 - x3_3*x5_0*k6_0 - x5_3*x3_0*k6_0 + x6_3*k4_0 + x6_3*k5_0 + x6_4,
		-x2_5 + 6838486782430686299035071132630075973697920250584089671969973910608713502903358704862830483072894667290471182617874180462193424855391356314223488979143494081738912668378445300235299223669454169180240866472469798832273931874771,
		-10*x1_3*x2_2*k1_0 - 10*x2_3*x1_2*k1_0 - 5*x1_4*x2_1*k1_0 - 5*x2_4*x1_1*k1_0 - x1_5*x2_0*k1_0 - x2_5*x1_0*k1_0 - x4_5*k2_0 - x4_5*k3_0 + x2_6,
		6*x2_2*x1_2*k1_0 + 4*x1_3*x2_1*k1_0 + 4*x2_3*x1_1*k1_0 + x1_4*x2_0*k1_0 + x2_4*x1_0*k1_0 - x4_4*k2_0 - x6_4*k4_0 + x1_5,
		-6*x2_2*x1_2*k1_0 - 4*x1_3*x2_1*k1_0 - 4*x2_3*x1_1*k1_0 - x1_4*x2_0*k1_0 - x2_4*x1_0*k1_0 + x4_4*k2_0 + x4_4*k3_0 + x4_5,
		-x3_5 - 9432538466279061969626251017997506200087012855382375310271219785484843806387873289962755838840610044593372993207503051477869979995636787726121110375505203077815042851356368062845396310043125746401457803519321295745845533921,
		10*x3_3*x5_2*k6_0 + 10*x5_3*x3_2*k6_0 + 5*x3_4*x5_1*k6_0 + 5*x5_4*x3_1*k6_0 + x3_5*x5_0*k6_0 + x5_5*x3_0*k6_0 - x4_5*k3_0 - x6_5*k5_0 + x3_6,
		6*x5_2*x3_2*k6_0 + 4*x3_3*x5_1*k6_0 + 4*x5_3*x3_1*k6_0 + x3_4*x5_0*k6_0 + x5_4*x3_0*k6_0 - x6_4*k4_0 - x6_4*k5_0 + x5_5,
		-6*x5_2*x3_2*k6_0 - 4*x3_3*x5_1*k6_0 - 4*x5_3*x3_1*k6_0 - x3_4*x5_0*k6_0 - x5_4*x3_0*k6_0 + x6_4*k4_0 + x6_4*k5_0 + x6_5,
		-x2_6 - 431823661101178583505815491042523134350206755479826028667592462756788009689941949869407738348139688220524420514575971797513660499210284796566738796263877621740295059112585851154193458024478996326631221756791028085954163425908665515752702603237782667975897109333515649,
		-x3_6 + 769870577908006120564448207008580735409610931139128690554050380585838323906530768958771900337367983660052503963075641299189963096098994417170599661737028858960886234625288502945817380331987105651194636010000131880991008488153208097961558067253161665808187776817089,
		z_aux - 1
    ]
end

# Critical pair degree sequence: 2,3,2,3,3,4,4,5,4,5,4,5,4,4,2,3,2,3,4,5,6,7
function Cholera(; np=AbstractAlgebra, internal_ordering=:degrevlex, k=np.QQ)
    R, (i_9,i_8,w_8,s_8,i_7,w_7,r_7,s_7,i_6,w_6,r_6,s_6,i_5,w_5,r_5,s_5,i_4,w_4,r_4,s_4,i_3,w_3,r_3,s_3,i_2,w_2,r_2,s_2,i_1,w_1,r_1,s_1,i_0,w_0,r_0,s_0,z_aux,mu_0,bi_0,bw_0,al_0,g_0,dz_0,k_0) = np.polynomial_ring(k, [:i_9,:i_8,:w_8,:s_8,:i_7,:w_7,:r_7,:s_7,:i_6,:w_6,:r_6,:s_6,:i_5,:w_5,:r_5,:s_5,:i_4,:w_4,:r_4,:s_4,:i_3,:w_3,:r_3,:s_3,:i_2,:w_2,:r_2,:s_2,:i_1,:w_1,:r_1,:s_1,:i_0,:w_0,:r_0,:s_0,:z_aux,:mu_0,:bi_0,:bw_0,:al_0,:g_0,:dz_0,:k_0], internal_ordering=internal_ordering)
    sys = [
    		-i_0 - r_0 - s_0 + 2980074650635745482722,
		r_0*mu_0 + r_0*al_0 - i_0*g_0 + r_1,
		-i_0*s_0*bi_0 - w_0*s_0*bw_0 + i_0*mu_0 + i_0*g_0 + i_1,
		i_0*s_0*bi_0 + w_0*s_0*bw_0 + s_0*mu_0 - r_0*al_0 + s_1 - mu_0,
		-i_0*k_0 + 187865766441801096442867884831449218259338,
		-i_1 - r_1 - s_1 - 587989940558935651497188412652797767056012,
		r_1*mu_0 + r_1*al_0 - i_1*g_0 + r_2,
		-s_1*i_0*bi_0 - i_1*s_0*bi_0 - s_1*w_0*bw_0 - w_1*s_0*bw_0 + i_1*mu_0 + i_1*g_0 + i_2,
		s_1*i_0*bi_0 + i_1*s_0*bi_0 + s_1*w_0*bw_0 + w_1*s_0*bw_0 + s_1*mu_0 - r_1*al_0 + s_2,
		-i_0*dz_0 + w_0*dz_0 + w_1,
		-i_1*k_0 + 1277410634018328126775017227746781373544979592248315878368595993560360406195112913062,
		-i_2*k_0 - 2341482585890381336538833296751745942818361578870820216616988459749244866760707031555541815643215948939122365114964260508063850,
		-2*i_1*s_1*bi_0 - s_2*i_0*bi_0 - i_2*s_0*bi_0 - 2*w_1*s_1*bw_0 - s_2*w_0*bw_0 - w_2*s_0*bw_0 + i_2*mu_0 + i_2*g_0 + i_3,
		-i_1*dz_0 + w_1*dz_0 + w_2,
		-i_3*k_0 + 3316015244139860530121857966267255563457696054341836040688382170569366288180772656573529715061868070249232375846729402189991897231243183810942121567633187898465260325146,
		-3*s_2*i_1*bi_0 - 3*i_2*s_1*bi_0 - s_3*i_0*bi_0 - i_3*s_0*bi_0 - 3*s_2*w_1*bw_0 - 3*w_2*s_1*bw_0 - s_3*w_0*bw_0 - w_3*s_0*bw_0 + i_3*mu_0 + i_3*g_0 + i_4,
		-i_2*dz_0 + w_2*dz_0 + w_3,
		2*i_1*s_1*bi_0 + s_2*i_0*bi_0 + i_2*s_0*bi_0 + 2*w_1*s_1*bw_0 + s_2*w_0*bw_0 + w_2*s_0*bw_0 + s_2*mu_0 - r_2*al_0 + s_3,
		-i_4*k_0 - 711767588783385145023388340846683054866008334206458998710644497577697279645891245332546655714533884848309778667396143632975186469663441981974068162334937266836639157223052701039453410718462860170791591818024586,
		-6*i_2*s_2*bi_0 - 4*s_3*i_1*bi_0 - 4*i_3*s_1*bi_0 - s_4*i_0*bi_0 - i_4*s_0*bi_0 - 6*w_2*s_2*bw_0 - 4*s_3*w_1*bw_0 - 4*w_3*s_1*bw_0 - s_4*w_0*bw_0 - w_4*s_0*bw_0 + i_4*mu_0 + i_4*g_0 + i_5,
		3*s_2*i_1*bi_0 + 3*i_2*s_1*bi_0 + s_3*i_0*bi_0 + i_3*s_0*bi_0 + 3*s_2*w_1*bw_0 + 3*w_2*s_1*bw_0 + s_3*w_0*bw_0 + w_3*s_0*bw_0 + s_3*mu_0 - r_3*al_0 + s_4,
		-i_3*dz_0 + w_3*dz_0 + w_4,
		r_2*mu_0 + r_2*al_0 - i_2*g_0 + r_3,
		-i_5*k_0 - 18665339724532000580604691285769607142920934357399551831480964903075741577870757797666897300030294582944931905607832352647880466738677411088168188229523985463421594620154643000945297999904581798490913330370177055753704211113126841752306417378998419358642,
		-10*s_3*i_2*bi_0 - 10*i_3*s_2*bi_0 - 5*s_4*i_1*bi_0 - 5*i_4*s_1*bi_0 - s_5*i_0*bi_0 - i_5*s_0*bi_0 - 10*s_3*w_2*bw_0 - 10*w_3*s_2*bw_0 - 5*s_4*w_1*bw_0 - 5*w_4*s_1*bw_0 - s_5*w_0*bw_0 - w_5*s_0*bw_0 + i_5*mu_0 + i_5*g_0 + i_6,
		6*i_2*s_2*bi_0 + 4*s_3*i_1*bi_0 + 4*i_3*s_1*bi_0 + s_4*i_0*bi_0 + i_4*s_0*bi_0 + 6*w_2*s_2*bw_0 + 4*s_3*w_1*bw_0 + 4*w_3*s_1*bw_0 + s_4*w_0*bw_0 + w_4*s_0*bw_0 + s_4*mu_0 - r_4*al_0 + s_5,
		-i_4*dz_0 + w_4*dz_0 + w_5,
		r_3*mu_0 + r_3*al_0 - i_3*g_0 + r_4,
		-i_6*k_0 + 83368010799974611909382386446631002045039777730921835098681336057320099362402100157422057780473244308796392262187575801000073554033141434465252932437061983204571177469369194997822796785158271078575363990002497784952965593543551975904451723114424207749293513848165974921635832137047353808469484746,
		-20*i_3*s_3*bi_0 - 15*s_4*i_2*bi_0 - 15*i_4*s_2*bi_0 - 6*s_5*i_1*bi_0 - 6*i_5*s_1*bi_0 - s_6*i_0*bi_0 - i_6*s_0*bi_0 - 20*w_3*s_3*bw_0 - 15*s_4*w_2*bw_0 - 15*w_4*s_2*bw_0 - 6*s_5*w_1*bw_0 - 6*w_5*s_1*bw_0 - s_6*w_0*bw_0 - w_6*s_0*bw_0 + i_6*mu_0 + i_6*g_0 + i_7,
		-i_5*dz_0 + w_5*dz_0 + w_6,
		10*s_3*i_2*bi_0 + 10*i_3*s_2*bi_0 + 5*s_4*i_1*bi_0 + 5*i_4*s_1*bi_0 + s_5*i_0*bi_0 + i_5*s_0*bi_0 + 10*s_3*w_2*bw_0 + 10*w_3*s_2*bw_0 + 5*s_4*w_1*bw_0 + 5*w_4*s_1*bw_0 + s_5*w_0*bw_0 + w_5*s_0*bw_0 + s_5*mu_0 - r_5*al_0 + s_6,
		r_4*mu_0 + r_4*al_0 - i_4*g_0 + r_5,
		-i_7*k_0 - 147967724216456888824593807280019292245302103850251837522159515303596838562643699163031711825578283990629163339503209909596209795604349100752179823072451003211338344645816712208045974936687218947527578753415691889254046093758180318902931670701150459085744978894957752687386282483147471744848096841479140419928757404958915652038052934379442,
		-35*s_4*i_3*bi_0 - 35*i_4*s_3*bi_0 - 21*s_5*i_2*bi_0 - 21*i_5*s_2*bi_0 - 7*s_6*i_1*bi_0 - 7*i_6*s_1*bi_0 - s_7*i_0*bi_0 - i_7*s_0*bi_0 - 35*s_4*w_3*bw_0 - 35*w_4*s_3*bw_0 - 21*s_5*w_2*bw_0 - 21*w_5*s_2*bw_0 - 7*s_6*w_1*bw_0 - 7*w_6*s_1*bw_0 - s_7*w_0*bw_0 - w_7*s_0*bw_0 + i_7*mu_0 + i_7*g_0 + i_8,
		-i_6*dz_0 + w_6*dz_0 + w_7,
		20*i_3*s_3*bi_0 + 15*s_4*i_2*bi_0 + 15*i_4*s_2*bi_0 + 6*s_5*i_1*bi_0 + 6*i_5*s_1*bi_0 + s_6*i_0*bi_0 + i_6*s_0*bi_0 + 20*w_3*s_3*bw_0 + 15*s_4*w_2*bw_0 + 15*w_4*s_2*bw_0 + 6*s_5*w_1*bw_0 + 6*w_5*s_1*bw_0 + s_6*w_0*bw_0 + w_6*s_0*bw_0 + s_6*mu_0 - r_6*al_0 + s_7,
		r_5*mu_0 + r_5*al_0 - i_5*g_0 + r_6,
		-i_8*k_0 - 674103273365030559391515752812072091829047312159829889448364059067220670276745456722507612194556136295436912111731469567590909041689687528599801557579800775473536986698752946179869493805440036196664100861484973318496238375134541202991386932229469343192215190091027687563318355352958089794174614458086269402784747884084002916914281085696398034501040776122991240187415034887078155914,
		-70*i_4*s_4*bi_0 - 56*s_5*i_3*bi_0 - 56*i_5*s_3*bi_0 - 28*s_6*i_2*bi_0 - 28*i_6*s_2*bi_0 - 8*s_7*i_1*bi_0 - 8*i_7*s_1*bi_0 - s_8*i_0*bi_0 - i_8*s_0*bi_0 - 70*w_4*s_4*bw_0 - 56*s_5*w_3*bw_0 - 56*w_5*s_3*bw_0 - 28*s_6*w_2*bw_0 - 28*w_6*s_2*bw_0 - 8*s_7*w_1*bw_0 - 8*w_7*s_1*bw_0 - s_8*w_0*bw_0 - w_8*s_0*bw_0 + i_8*mu_0 + i_8*g_0 + i_9,
		35*s_4*i_3*bi_0 + 35*i_4*s_3*bi_0 + 21*s_5*i_2*bi_0 + 21*i_5*s_2*bi_0 + 7*s_6*i_1*bi_0 + 7*i_6*s_1*bi_0 + s_7*i_0*bi_0 + i_7*s_0*bi_0 + 35*s_4*w_3*bw_0 + 35*w_4*s_3*bw_0 + 21*s_5*w_2*bw_0 + 21*w_5*s_2*bw_0 + 7*s_6*w_1*bw_0 + 7*w_6*s_1*bw_0 + s_7*w_0*bw_0 + w_7*s_0*bw_0 + s_7*mu_0 - r_7*al_0 + s_8,
		-i_7*dz_0 + w_7*dz_0 + w_8,
		r_6*mu_0 + r_6*al_0 - i_6*g_0 + r_7,
		-i_2 - r_2 - s_2 + 116014600548595291886176008075760517732197215875621423410450064,
		-i_3 - r_3 - s_3 - 22890506473045791995115165087583625844781809884223370988178227721730602360097763008,
		-i_4 - r_4 - s_4 + 4516459860352427093382687739390813045476029732369311452081971342192413898927916571756156200650811781376,
		-i_5 - r_5 - s_5 - 891129678331685668505471310021741105387507255043152430038674460356075455664492139206771807662640426076290683418128023964672,
		-i_6 - r_6 - s_6 + 175826228541211399892699778334585423146889420531107175406069200855246858982303884279890012859457287853465355668537083211173082437423411127619584,
		-i_7 - r_7 - s_7 - 34691766411486904318482481905083491952664034705323971165028302785858442473221550597346325692604524560056453654749755021785868967493673755968388313225830098258280448,
		-i_9*k_0 + 7471136511829128504299483090560499068662192052395856667723160028978279590694775782794424094143010609362277289590643587793496605382946287380413738801546605018611179544663229484501796360121821537396946476224658214561900384576832669227906100996542315194986195519446703769836611908531577602122585984958790071117884682479803541091090718034761019332672999051148441055356488226926951744455387486628157925225943716368750244978261890,
		z_aux - 1
    ]
end

# Critical pair degree sequence: 4,5,5,6,7,6,7,8,7,8,9,8,9,9,9,10,7,8,9,10,9,10,9,9,10,10,10,9,10,9,9,10,9,9,8,7,4,4,3,4,2,3,4,5,4,6,4,5,6,4,5,6,7,6,7,5,8,9,10,11,12,13,14,15,16,17,18,19
function Goodwin_with_weights(; np=AbstractAlgebra, internal_ordering=:degrevlex, k=np.QQ)
    R, (x1_9,x1_8,x4_8,x1_7,x2_7,x3_7,x4_7,x1_6,x2_6,x3_6,x4_6,x1_5,x2_5,x3_5,x4_5,x1_4,x2_4,x3_4,x4_4,x1_3,x2_3,x3_3,x4_3,x1_2,x2_2,x3_2,x4_2,x1_1,x2_1,x3_1,x4_1,x1_0,x2_0,x3_0,x4_0,z_aux,b_0,c_0,alpha_0,beta_0,gama_0,delta_0,sigma_0) = np.polynomial_ring(k, [:x1_9,:x1_8,:x4_8,:x1_7,:x2_7,:x3_7,:x4_7,:x1_6,:x2_6,:x3_6,:x4_6,:x1_5,:x2_5,:x3_5,:x4_5,:x1_4,:x2_4,:x3_4,:x4_4,:x1_3,:x2_3,:x3_3,:x4_3,:x1_2,:x2_2,:x3_2,:x4_2,:x1_1,:x2_1,:x3_1,:x4_1,:x1_0,:x2_0,:x3_0,:x4_0,:z_aux,:b_0,:c_0,:alpha_0,:beta_0,:gama_0,:delta_0,:sigma_0], internal_ordering=internal_ordering)
    sys = [
    		-x1_0 + 66626618109538133133,
		x1_0*x4_0^2*b_0 + x1_1*x4_0^2 + x1_0*b_0*c_0 + x1_1*c_0 - 1,
		-x1_1 - k(109365372900331087094671182948099025791485844124302772199831)//k(20532014486747450052),
		x4_1^2*x1_0*b_0 + x1_1*x4_0^2*b_0 + x1_1*x4_1^2 + x1_2*x4_0^2 + x1_1*b_0*c_0 + x1_2*c_0,
		-x2_0^3*x4_0^2*gama_0*sigma_0 + x3_0^3*x4_0^2*delta_0*sigma_0 + x4_1^2*x3_0^3,
		-x1_2 + k(700046698782145700339511765019070604588317901871928458981218413415494465333157989305079273556783733467711209787417995)//k(1643910726782351006692742235907081361548497842545856874048),
		2*x1_1*x4_1^2*b_0 + x4_2^2*x1_0*b_0 + x1_2*x4_0^2*b_0 + x4_2^2*x1_1 + 2*x1_2*x4_1^2 + x1_3*x4_0^2 + x1_2*b_0*c_0 + x1_3*c_0,
		-x4_1^2*x2_0^3*gama_0*sigma_0 - x2_1^3*x4_0^2*gama_0*sigma_0 + x4_1^2*x3_0^3*delta_0*sigma_0 + x3_1^3*x4_0^2*delta_0*sigma_0 + x3_1^3*x4_1^2 + x4_2^2*x3_0^3,
		x2_0^3*beta_0^4 + x2_1^3 - x1_0*alpha_0,
		-x2_0^3*gama_0 + x3_0^3*delta_0 + x3_1^3,
		-x1_3 - k(746832029613852051468613943547303704721851393864042529307888871434982983597843743482981840300040737111171058150705189565199900965660717911805427355150537058230893366173423925)//k(21936818712831695217261846127345386652907426055379126084735493887806963489926767660633850107392),
		3*x4_2^2*x1_1*b_0 + 3*x1_2*x4_1^2*b_0 + x4_3^2*x1_0*b_0 + x1_3*x4_0^2*b_0 + 3*x1_2*x4_2^2 + x4_3^2*x1_1 + 3*x1_3*x4_1^2 + x1_4*x4_0^2 + x1_3*b_0*c_0 + x1_4*c_0,
		-2*x2_1^3*x4_1^2*gama_0*sigma_0 - x4_2^2*x2_0^3*gama_0*sigma_0 - x2_2^3*x4_0^2*gama_0*sigma_0 + 2*x3_1^3*x4_1^2*delta_0*sigma_0 + x4_2^2*x3_0^3*delta_0*sigma_0 + x3_2^3*x4_0^2*delta_0*sigma_0 + 2*x4_2^2*x3_1^3 + x3_2^3*x4_1^2 + x4_3^2*x3_0^3,
		x2_1^3*beta_0^4 + x2_2^3 - x1_1*alpha_0,
		-x2_1^3*gama_0 + x3_1^3*delta_0 + x3_2^3,
		-x1_4 + k(14400496282353147779240440064051517594542619061482302119481500918335678783139767285489155810529354387115922947240349381533886315899816080905890011513422804140647502909509376093755145644195474010268903020157172111149233576536203218793)//k(5269162207650950409816206966484017258892742291065923526300892532467949058939610946177006000036951932644775528967156598985881002778624),
		6*x1_2*x4_2^2*b_0 + 4*x4_3^2*x1_1*b_0 + 4*x1_3*x4_1^2*b_0 + x4_4^2*x1_0*b_0 + x1_4*x4_0^2*b_0 + 4*x4_3^2*x1_2 + 6*x1_3*x4_2^2 + x4_4^2*x1_1 + 4*x1_4*x4_1^2 + x1_5*x4_0^2 + x1_4*b_0*c_0 + x1_5*c_0,
		-3*x4_2^2*x2_1^3*gama_0*sigma_0 - 3*x2_2^3*x4_1^2*gama_0*sigma_0 - x4_3^2*x2_0^3*gama_0*sigma_0 - x2_3^3*x4_0^2*gama_0*sigma_0 + 3*x4_2^2*x3_1^3*delta_0*sigma_0 + 3*x3_2^3*x4_1^2*delta_0*sigma_0 + x4_3^2*x3_0^3*delta_0*sigma_0 + x3_3^3*x4_0^2*delta_0*sigma_0 + 3*x3_2^3*x4_2^2 + 3*x4_3^2*x3_1^3 + x3_3^3*x4_1^2 + x4_4^2*x3_0^3,
		-x2_2^3*gama_0 + x3_2^3*delta_0 + x3_3^3,
		x2_2^3*beta_0^4 + x2_3^3 - x1_2*alpha_0,
		-x1_5 + k(6023146484867968982831395072975853459840428675793078655027031200915625335251026264995591788314817061012960143127465896778104301317190807787359147785084035892219997768561722127990674883087817275061642343495909207835010250562445419129856000180674740973000608844919335360292523100181672948353464614971220671779)//k(105469829566706030449828814232727792805929839103607170061003122354503588455876618384576557910763445930104623404066672449815409164413677718859087964596894741231747453190144),
		10*x4_3^2*x1_2*b_0 + 10*x1_3*x4_2^2*b_0 + 5*x4_4^2*x1_1*b_0 + 5*x1_4*x4_1^2*b_0 + x4_5^2*x1_0*b_0 + x1_5*x4_0^2*b_0 + 10*x1_3*x4_3^2 + 5*x4_4^2*x1_2 + 10*x1_4*x4_2^2 + x4_5^2*x1_1 + 5*x1_5*x4_1^2 + x1_6*x4_0^2 + x1_5*b_0*c_0 + x1_6*c_0,
		-6*x2_2^3*x4_2^2*gama_0*sigma_0 - 4*x4_3^2*x2_1^3*gama_0*sigma_0 - 4*x2_3^3*x4_1^2*gama_0*sigma_0 - x4_4^2*x2_0^3*gama_0*sigma_0 - x2_4^3*x4_0^2*gama_0*sigma_0 + 6*x3_2^3*x4_2^2*delta_0*sigma_0 + 4*x4_3^2*x3_1^3*delta_0*sigma_0 + 4*x3_3^3*x4_1^2*delta_0*sigma_0 + x4_4^2*x3_0^3*delta_0*sigma_0 + x3_4^3*x4_0^2*delta_0*sigma_0 + 6*x4_3^2*x3_2^3 + 4*x3_3^3*x4_2^2 + 4*x4_4^2*x3_1^3 + x3_4^3*x4_1^2 + x4_5^2*x3_0^3,
		-x2_3^3*gama_0 + x3_3^3*delta_0 + x3_4^3,
		x2_3^3*beta_0^4 + x2_4^3 - x1_3*alpha_0,
		-x1_6 + k(253949913262840323088574863777035021031954734282314408183787308682017877732179958729999565627829424377744430379593102942481557434380050527427389510359095048100556631522108698594500953750237585970360542368910419442925329493834675973473913093385116882696735396338749992298879460076980869416976912998032080964055940429022340750367620324990842004351034250664637608520392629227630991313557)//k(8444518889684488176011846605598492563045475803377931581512194493985531879976633361612792818682212573986093149138098413137707835069293281368636953083499711979427409877779945874471845453307638156373674078765056),
		20*x1_3*x4_3^2*b_0 + 15*x4_4^2*x1_2*b_0 + 15*x1_4*x4_2^2*b_0 + 6*x4_5^2*x1_1*b_0 + 6*x1_5*x4_1^2*b_0 + x4_6^2*x1_0*b_0 + x1_6*x4_0^2*b_0 + 15*x4_4^2*x1_3 + 20*x1_4*x4_3^2 + 6*x4_5^2*x1_2 + 15*x1_5*x4_2^2 + x4_6^2*x1_1 + 6*x1_6*x4_1^2 + x1_7*x4_0^2 + x1_6*b_0*c_0 + x1_7*c_0,
		-10*x4_3^2*x2_2^3*gama_0*sigma_0 - 10*x2_3^3*x4_2^2*gama_0*sigma_0 - 5*x4_4^2*x2_1^3*gama_0*sigma_0 - 5*x2_4^3*x4_1^2*gama_0*sigma_0 - x4_5^2*x2_0^3*gama_0*sigma_0 - x2_5^3*x4_0^2*gama_0*sigma_0 + 10*x4_3^2*x3_2^3*delta_0*sigma_0 + 10*x3_3^3*x4_2^2*delta_0*sigma_0 + 5*x4_4^2*x3_1^3*delta_0*sigma_0 + 5*x3_4^3*x4_1^2*delta_0*sigma_0 + x4_5^2*x3_0^3*delta_0*sigma_0 + x3_5^3*x4_0^2*delta_0*sigma_0 + 10*x3_3^3*x4_3^2 + 10*x4_4^2*x3_2^3 + 5*x3_4^3*x4_2^2 + 5*x4_5^2*x3_1^3 + x3_5^3*x4_1^2 + x4_6^2*x3_0^3,
		x2_4^3*beta_0^4 + x2_5^3 - x1_4*alpha_0,
		-x2_4^3*gama_0 + x3_4^3*delta_0 + x3_5^3,
		-x1_7 - k(121666364780579037201618891345504960968107797254699155205895968421189028955557603093854243603214720972195545210999472311147263143262522705177829447965841688225081819816791643225692481042571672286399621474800225009939636012375612436733407971481545257067325917538109866979012665403604268770054516144410900784075735332150384593144957503364412194553880477271088774944909630657666985161908907917526367546169334191425080570985591466738945720763238706673089596743029443)//k(338058284398464342634237460051035882962821424940867702737940688189897411835746130255995754001942685562871984955545276090349686810895527222299247607721903370611530697697112014593813084860842123183476550820848674364930866131300064595142068112719872),
		35*x4_4^2*x1_3*b_0 + 35*x1_4*x4_3^2*b_0 + 21*x4_5^2*x1_2*b_0 + 21*x1_5*x4_2^2*b_0 + 7*x4_6^2*x1_1*b_0 + 7*x1_6*x4_1^2*b_0 + x4_7^2*x1_0*b_0 + x1_7*x4_0^2*b_0 + 35*x1_4*x4_4^2 + 21*x4_5^2*x1_3 + 35*x1_5*x4_3^2 + 7*x4_6^2*x1_2 + 21*x1_6*x4_2^2 + x4_7^2*x1_1 + 7*x1_7*x4_1^2 + x1_8*x4_0^2 + x1_7*b_0*c_0 + x1_8*c_0,
		-20*x2_3^3*x4_3^2*gama_0*sigma_0 - 15*x4_4^2*x2_2^3*gama_0*sigma_0 - 15*x2_4^3*x4_2^2*gama_0*sigma_0 - 6*x4_5^2*x2_1^3*gama_0*sigma_0 - 6*x2_5^3*x4_1^2*gama_0*sigma_0 - x4_6^2*x2_0^3*gama_0*sigma_0 - x2_6^3*x4_0^2*gama_0*sigma_0 + 20*x3_3^3*x4_3^2*delta_0*sigma_0 + 15*x4_4^2*x3_2^3*delta_0*sigma_0 + 15*x3_4^3*x4_2^2*delta_0*sigma_0 + 6*x4_5^2*x3_1^3*delta_0*sigma_0 + 6*x3_5^3*x4_1^2*delta_0*sigma_0 + x4_6^2*x3_0^3*delta_0*sigma_0 + x3_6^3*x4_0^2*delta_0*sigma_0 + 20*x4_4^2*x3_3^3 + 15*x3_4^3*x4_3^2 + 15*x4_5^2*x3_2^3 + 6*x3_5^3*x4_2^2 + 6*x4_6^2*x3_1^3 + x3_6^3*x4_1^2 + x4_7^2*x3_0^3,
		-x2_5^3*gama_0 + x3_5^3*delta_0 + x3_6^3,
		x2_5^3*beta_0^4 + x2_6^3 - x1_5*alpha_0,
		-x1_8 - k(1035478426339532804280184425675865212583192915679354737011282251768492971659500295192883199178837870042316593767404434135839502034165360094555662944021415613923971520828829759859200678256806858872793173201520764955351575794871599304992700879752149237346260831272097084348101519836269315550398293457515278349705693029148551716663479630565103093831227504859508525747907421333476785031504646003998709330663881303776233995302032290120641287277558916223602660599661879631991614501351459525563588307585853179639238170094472370248825687845888525)//k(1002477159211585928057936315784867313939782158331411691253395651104546961303476537643951251274759542426673822147953537622807982818892001485731340326345370249669999250300655893874635217205774850859399132984297060844138082946659677229672722667405027229797464480024009081242726447448064),
		70*x1_4*x4_4^2*b_0 + 56*x4_5^2*x1_3*b_0 + 56*x1_5*x4_3^2*b_0 + 28*x4_6^2*x1_2*b_0 + 28*x1_6*x4_2^2*b_0 + 8*x4_7^2*x1_1*b_0 + 8*x1_7*x4_1^2*b_0 + x4_8^2*x1_0*b_0 + x1_8*x4_0^2*b_0 + 56*x4_5^2*x1_4 + 70*x1_5*x4_4^2 + 28*x4_6^2*x1_3 + 56*x1_6*x4_3^2 + 8*x4_7^2*x1_2 + 28*x1_7*x4_2^2 + x4_8^2*x1_1 + 8*x1_8*x4_1^2 + x1_9*x4_0^2 + x1_8*b_0*c_0 + x1_9*c_0,
		-35*x4_4^2*x2_3^3*gama_0*sigma_0 - 35*x2_4^3*x4_3^2*gama_0*sigma_0 - 21*x4_5^2*x2_2^3*gama_0*sigma_0 - 21*x2_5^3*x4_2^2*gama_0*sigma_0 - 7*x4_6^2*x2_1^3*gama_0*sigma_0 - 7*x2_6^3*x4_1^2*gama_0*sigma_0 - x4_7^2*x2_0^3*gama_0*sigma_0 - x2_7^3*x4_0^2*gama_0*sigma_0 + 35*x4_4^2*x3_3^3*delta_0*sigma_0 + 35*x3_4^3*x4_3^2*delta_0*sigma_0 + 21*x4_5^2*x3_2^3*delta_0*sigma_0 + 21*x3_5^3*x4_2^2*delta_0*sigma_0 + 7*x4_6^2*x3_1^3*delta_0*sigma_0 + 7*x3_6^3*x4_1^2*delta_0*sigma_0 + x4_7^2*x3_0^3*delta_0*sigma_0 + x3_7^3*x4_0^2*delta_0*sigma_0 + 35*x3_4^3*x4_4^2 + 35*x4_5^2*x3_3^3 + 21*x3_5^3*x4_3^2 + 21*x4_6^2*x3_2^3 + 7*x3_6^3*x4_2^2 + 7*x4_7^2*x3_1^3 + x3_7^3*x4_1^2 + x4_8^2*x3_0^3,
		-x2_6^3*gama_0 + x3_6^3*delta_0 + x3_7^3,
		x2_6^3*beta_0^4 + x2_7^3 - x1_6*alpha_0,
		-x1_9 + k(74498781299447588606588752974911930670277856803642844339909705301308548240554872329687775913191115501078578240897919561056135671823288981781246290384559929160172441876921532712805754611079471528573837636495175797816357679171877191693241936172759683737163552555041129791142697756563432473275798260805176836905469839505534348763789895151148761693528397708899825983488347397557503480751294060701474762176019247748646004831690933150503768771039290159708433683307227220068724618779227639372844273945676326911995918286827612406310569470979586395325145160102199899416795812552081561134366302548834799830328753962802162295)//k(30099024558317116707323832373444070201560469097079658125143427669949259251197216121216796228195086343457976491681632586935429455336434256495349995437638886080316189055707152086758344142451695886597267848791450868314255912444883882862663815068204672969810371079224322972046994652220833761909795444682289947002438393266176),
		x3_0^3*x4_0^2*z_aux + x3_0^3*z_aux*c_0 - 1
    ]
end

###
# Random generation

# A random polynomial with integer coefficients.
function random_poly(rng::AbstractRNG, ring, maxdeg::Int, nterms::Int, coeffsz::Integer)
    rand(rng, ring, 0:maxdeg, 1:nterms, 1:coeffsz)
end

# A random polynomial system over k.
function random_generating_set(
        rng::AbstractRNG,
        k,
        ord::Symbol,
        nvars::Int,
        maxdeg::Int,
        nterms::Int,
        npolys::Int,
        coeffsz::Integer;
        np=AbstractAlgebra
    )
    R, _ = np.polynomial_ring(np.ZZ, ["x$i" for i in 1:nvars], internal_ordering=ord)
    polys = [random_poly(rng, R, maxdeg, nterms, coeffsz) for _ in 1:npolys]
    polys = filter(!iszero, polys)
    isempty(polys) && return polys
    polys = map(f -> np.map_coefficients(c -> k(c), f), polys)
    polys
end

end # module Examples

using .Examples
