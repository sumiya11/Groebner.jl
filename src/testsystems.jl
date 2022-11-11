# Systems used for testing and benchmarking

# The cyclic-n system
# (not to be confused with root-n system)!!
#
# Göran Björck and Ralf Fröberg: 
# `A faster way to count the solutions of inhomogeneous systems of algebraic equations, with applications to cyclic n-roots', 
# in J. Symbolic Computation (1991) 12, pp 329–336.
function cyclicn(n; np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, z = np.PolynomialRing(ground, ["z$i" for i in 1:n], ordering=ordering)
    [
        (sum(prod(z[(k-1) % n + 1] for k in j:j+m) for j in 1:n) 
        for m=0:(n-2))...,prod(z)-1
    ]
end

#------------------------------------------------------------------------------

# The root-n system
# (not to be confused with cyclic-n system)!!
function rootn(n; np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, xs = np.PolynomialRing(
        ground, 
        ["x$i" for i in 1:n], 
        ordering=ordering
    )
    ans = [
        sum(map(prod, Combinatorics.combinations(xs, i)))
        for i in 1:n
    ]
    ans[end] -= (-1)^(n - 1)
    ans
end

#------------------------------------------------------------------------------

# The reimer-n system
function reimern(n; np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, xs = np.PolynomialRing(ground, ["x$i" for i in 1:n])
    [
        sum((-1)^(i+1)*2*xs[i]^j for i in 1:n) - 1
        for j in 2:(n+1)
    ]
end

#------------------------------------------------------------------------------

# The katsura-n system
#
# S. Katsura, W. Fukuda, S. Inawashiro, N.M. Fujiki and R. Gebauer,
#  Cell Biophysics, Vol 11, pages 309–319, 1987.
function katsuran(n; np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, x = np.PolynomialRing(ground, ["x$i" for i in 0:n], ordering=ordering)
    [
        (sum(x[abs(l)+1]*x[abs(m-l)+1] for l=-n:n if abs(m-l)<=n) -
        x[m+1] for m=0:n-1)...,
        x[1] + 2sum(x[i+1] for i=1:n) - 1
    ]
end

#------------------------------------------------------------------------------

# The noon-n system
# 
function noonn(n; np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    without(x, k) = x[1:end .!= k]

    R, xs = np.PolynomialRing(ground, ["x$i" for i in 1:n], ordering=ordering)
    fs = zeros(R, n)
    for i in 1:n
        other = without(xs, i)
        fs[i] = xs[i] * (10*sum(other .^ 2) - 11) + 10
    end
    fs
end

#------------------------------------------------------------------------------

# The henrion-5 system
#
# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion5.ms
function henrion5(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (f1,f2,f3,f4,f5,t) = np.PolynomialRing(ground, ["f1","f2","f3","f4","f5","t"], ordering=ordering)
    [
        2*f1*f2*f3*f4*f5-9823275,
        ground(21)//5*f1*f2*f4*f5+ground(16)//5*f1*f3*f4*f5+ground(9)//5*f2*f3*f4*f5+ground(24)//5*f1*f2*f3*f5+5*f4*f3*f1*f2-4465125,
        ground(14)//5*f4*f5*f1+ground(14)//5*f4*f5*f2+ground(8)//5*f3*f4*f5+ground(18)//5*f1*f2*f5+ground(24)//5*f1*f3*f5+ground(18)//5*f2*f3*f5+4*f3*f1*f2+6*f1*f2*f4+6*f3*f4*f1+4*f2*f3*f4-441486,
        ground(7)//5*f4*f5+ground(12)//5*f5*f1+ground(12)//5*f5*f2+ground(12)//5*f5*f3+3*f1*f2+4*f3*f1+4*f4*f1+3*f2*f3+4*f4*f2+3*f3*f4-15498,
        ground(6)//5*f5+2*f4+2*f3+2*f2+2*f1-215,
        f1+2*f2+3*f3+4*f4+5*f5+6*t
    ]
end

# The henrion-6 system
#
# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion6.ms
function henrion6(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (f1,f2,f3,f4,f5,f6) = np.PolynomialRing(ground, ["f1","f2","f3","f4","f5","f6"], ordering=ordering)

    [
        2*f1*f2*f3*f4*f5*f6-1404728325,
        6*f5*f4*f3*f1*f2+ground(11)//6*f2*f3*f4*f5*f6+ground(16)//3*f1*f2*f3*f5*f6+ground(9)//2*f1*f2*f4*f5*f6+ground(10)//3*f1*f3*f4*f5*f6+ground(35)//6*f1*f2*f3*f4*f6-648336150,
        5*f4*f3*f1*f2+5*f2*f3*f4*f5+ground(5)//3*f3*f4*f5*f6+8*f1*f2*f3*f5+9*f1*f2*f4*f5+8*f1*f3*f4*f5+4*f1*f2*f5*f6+ground(16)//3*f1*f3*f5*f6+3*f1*f4*f5*f6+4*f2*f3*f5*f6+3*f2*f4*f5*f6+ground(14)//3*f1*f2*f3*f6+7*f1*f2*f4*f6+7*f1*f3*f4*f6+ground(14)//3*f2*f3*f4*f6-67597623,
        6*f1*f2*f5+8*f1*f3*f5+6*f2*f3*f5+ground(8)//3*f5*f6*f3+ground(8)//3*f5*f6*f2+ground(8)//3*f5*f6*f1+ground(7)//2*f1*f2*f6+ground(14)//3*f1*f3*f6+ground(14)//3*f1*f4*f6+ground(7)//2*f2*f3*f6+ground(14)//3*f2*f4*f6+ground(7)//2*f3*f4*f6+6*f4*f5*f1+ground(3)//2*f4*f5*f6+4*f3*f1*f2+4*f2*f3*f4+6*f3*f4*f1+4*f3*f4*f5+6*f1*f2*f4+6*f4*f5*f2-2657700,
        ground(4)//3*f5*f6+ground(7)//3*f6*f1+ground(7)//3*f6*f2+ground(7)//3*f6*f3+ground(7)//3*f6*f4+3*f1*f2+4*f3*f1+4*f4*f1+4*f5*f1+3*f2*f3+4*f4*f2+4*f5*f2+3*f3*f4+4*f5*f3+3*f4*f5-46243,
        ground(7)//6*f6+2*f5+2*f4+2*f3+2*f2+2*f1-358
    ]
end

# The henrion-7 system
#
# Source:
# https://gitlab.lip6.fr/eder/msolve-examples/-/raw/master/zero-dimensional/henrion7.ms
function henrion7(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (f1,f2,f3,f4,f5,f6,f7) = np.PolynomialRing(ground, ["f1","f2","f3","f4","f5","f6","f7"], ordering=ordering)
    [
        2*f1*f2*f3*f4*f5*f6*f7-273922023375,
        ground(45)//7*f1*f2*f3*f4*f6*f7+ground(40)//7*f1*f2*f3*f5*f6*f7+ground(33)//7*f1*f2*f4*f5*f6*f7+ground(24)//7*f1*f3*f4*f5*f6*f7+ground(48)//7*f1*f2*f3*f4*f5*f7+7*f6*f5*f4*f3*f1*f2+ground(13)//7*f2*f3*f4*f5*f6*f7-127830277575,
        6*f5*f4*f3*f1*f2+6*f2*f3*f4*f5*f6+ground(12)//7*f3*f4*f5*f6*f7+10*f1*f2*f3*f4*f6+12*f1*f2*f3*f5*f6+12*f1*f2*f4*f5*f6+10*f1*f3*f4*f5*f6+ground(36)//7*f1*f2*f3*f6*f7+ground(54)//7*f1*f2*f4*f6*f7+ground(30)//7*f1*f2*f5*f6*f7+ground(54)//7*f1*f3*f4*f6*f7+ground(40)//7*f1*f3*f5*f6*f7+ground(22)//7*f1*f4*f5*f6*f7+ground(36)//7*f2*f3*f4*f6*f7+ground(30)//7*f2*f3*f5*f6*f7+ground(22)//7*f2*f4*f5*f6*f7+ground(40)//7*f1*f2*f3*f4*f7+ground(64)//7*f1*f2*f3*f5*f7+ground(72)//7*f1*f2*f4*f5*f7+ground(64)//7*f1*f3*f4*f5*f7+ground(40)//7*f2*f3*f4*f5*f7-13829872635,
        -585849123+5*f4*f3*f1*f2+5*f2*f3*f4*f5+5*f3*f4*f5*f6+8*f1*f2*f3*f5+9*f1*f2*f4*f5+8*f1*f3*f4*f5+ground(11)//7*f4*f5*f6*f7+8*f1*f2*f3*f6+12*f1*f2*f4*f6+9*f1*f2*f5*f6+12*f1*f3*f4*f6+12*f1*f3*f5*f6+8*f1*f4*f5*f6+8*f2*f3*f4*f6+9*f2*f3*f5*f6+8*f2*f4*f5*f6+ground(27)//7*f1*f2*f6*f7+ground(36)//7*f1*f3*f6*f7+ground(36)//7*f1*f4*f6*f7+ground(20)//7*f1*f5*f6*f7+ground(27)//7*f2*f3*f6*f7+ground(36)//7*f2*f4*f6*f7+ground(20)//7*f2*f5*f6*f7+ground(27)//7*f3*f4*f6*f7+ground(20)//7*f3*f5*f6*f7+ground(32)//7*f1*f2*f3*f7+ground(48)//7*f1*f2*f4*f7+ground(48)//7*f1*f2*f5*f7+ground(48)//7*f1*f3*f4*f7+ground(64)//7*f1*f3*f5*f7+ground(48)//7*f1*f4*f5*f7+ground(32)//7*f2*f3*f4*f7+ground(48)//7*f2*f3*f5*f7+ground(48)//7*f2*f4*f5*f7+ground(32)//7*f3*f4*f5*f7,
        -11675085+6*f1*f2*f6+8*f1*f3*f6+8*f1*f4*f6+6*f2*f3*f6+8*f2*f4*f6+6*f3*f4*f6+ground(18)//7*f6*f7*f4+ground(18)//7*f6*f7*f3+ground(18)//7*f6*f7*f2+ground(18)//7*f6*f7*f1+ground(24)//7*f1*f2*f7+ground(32)//7*f1*f3*f7+ground(32)//7*f1*f4*f7+ground(32)//7*f1*f5*f7+ground(24)//7*f2*f3*f7+ground(32)//7*f2*f4*f7+ground(32)//7*f2*f5*f7+ground(24)//7*f3*f4*f7+ground(32)//7*f3*f5*f7+ground(24)//7*f4*f5*f7+
        ground(10)//7*f5*f6*f7+6*f1*f2*f5+8*f1*f3*f5+6*f2*f3*f5+6*f5*f6*f3+6*f5*f6*f2+6*f5*f6*f1+6*f1*f2*f4+6*f4*f5*f2+6*f4*f5*f1+4*f4*f5*f6+6*f3*f4*f1+4*f3*f4*f5+4*f3*f1*f2+4*f2*f3*f4,
        ground(9)//7*f6*f7+ground(16)//7*f7*f1+ground(16)//7*f7*f2+ground(16)//7*f7*f3+ground(16)//7*f7*f4+ground(16)//7*f7*f5+3*f1*f2+4*f3*f1+4*f4*f1+4*f5*f1+4*f6*f1+3*f2*f3+4*f4*f2+4*f5*f2+4*f6*f2+3*f3*f4+4*f5*f3+4*f6*f3+3*f4*f5+4*f6*f4+3*f5*f6-116053,
        ground(8)//7*f7+2*f6+2*f5+2*f4+2*f3+2*f2+2*f1-553
    ]
end

#------------------------------------------------------------------------------

# Eco-n systems:

function eco5(;np=AbstractAlgebra, ground=np.GF(2^31-1), ordering=:lex)
    _, (x1, x2, x3, x4, x5) = np.PolynomialRing(ground, ["x$i" for i in 1:5], ordering=ordering)
    [
    (x1 + x1*x2 + x2*x3 + x3*x4)*x5 - 1,
     (x2 + x1*x3 + x2*x4)*x5 - 2,
             (x3 + x1*x4)*x5 - 3,
                       x4*x5 - 4,
           x1 + x2 + x3 + x4 + 1
    ]
end

function eco7(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x1, x2, x3, x4, x5, x6, x7) = np.PolynomialRing(ground, ["x$i" for i in 1:7], ordering=ordering)
    [
        (x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6)*x7 - 1,
     (x2 + x1*x3 + x2*x4 + x3*x5 + x4*x6)*x7 - 2,
     (x3 + x1*x4 + x2*x5 + x3*x6)*x7 - 3,
     (x4 + x1*x5 + x2*x6)*x7 - 4,
     (x5 + x1*x6)*x7 - 5,
     x6*x7 - 6,
     x1 + x2 + x3 + x4 + x5 + x6 + 1
    ]
end

function eco10(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9) = np.PolynomialRing(ground, ["x$i" for i in 1:11], ordering=ordering)
    [
        x0*x1*x9+x1*x2*x9+x2*x3*x9+x3*x4*x9+x4*x5*x9+x5*x6*x9+x6*x7*x9+x7*x8*x9+x0*x9-1,
        x0*x2*x9+x1*x3*x9+x2*x4*x9+x3*x5*x9+x4*x6*x9+x5*x7*x9+x6*x8*x9+x1*x9-2,
        x0*x3*x9+x1*x4*x9+x2*x5*x9+x3*x6*x9+x4*x7*x9+x5*x8*x9+x2*x9-3,
        x0*x4*x9+x1*x5*x9+x2*x6*x9+x3*x7*x9+x4*x8*x9+x3*x9-4,
        x0*x5*x9+x1*x6*x9+x2*x7*x9+x3*x8*x9+x4*x9-5,
        x0*x6*x9+x1*x7*x9+x2*x8*x9+x5*x9-6,
        x0*x7*x9+x1*x8*x9+x6*x9-7,
        x0*x8*x9+x7*x9-8,
        x8*x9-9,
        x0+x1+x2+x3+x4+x5+x6+x7+x8+1
    ]
end

function eco11(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = np.PolynomialRing(ground, ["x$i" for i in 1:11], ordering=ordering)
    [
        x0*x1*x10+x1*x2*x10+x2*x3*x10+x3*x4*x10+x4*x5*x10+x5*x6*x10+x6*x7*x10+x7*x8*x10+x8*x9*x10+x0*x10-1,
        x0*x2*x10+x1*x3*x10+x2*x4*x10+x3*x5*x10+x4*x6*x10+x5*x7*x10+x6*x8*x10+x7*x9*x10+x1*x10-2,
        x0*x3*x10+x1*x4*x10+x2*x5*x10+x3*x6*x10+x4*x7*x10+x5*x8*x10+x6*x9*x10+x2*x10-3,
        x0*x4*x10+x1*x5*x10+x2*x6*x10+x3*x7*x10+x4*x8*x10+x5*x9*x10+x3*x10-4,
        x0*x5*x10+x1*x6*x10+x2*x7*x10+x3*x8*x10+x4*x9*x10+x4*x10-5,
        x0*x6*x10+x1*x7*x10+x2*x8*x10+x3*x9*x10+x5*x10-6,
        x0*x7*x10+x1*x8*x10+x2*x9*x10+x6*x10-7,
        x0*x8*x10+x1*x9*x10+x7*x10-8,
        x0*x9*x10+x8*x10-9,
        x9*x10-10,
        x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+1
    ]
end

function eco12(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) = np.PolynomialRing(ground, ["x$i" for i in 1:12], ordering=ordering)
    [
        x0*x1*x11+x1*x2*x11+x2*x3*x11+x3*x4*x11+x4*x5*x11+x5*x6*x11+x6*x7*x11+x7*x8*x11+x8*x9*x11+x9*x10*x11+x0*x11-1,
        x0*x2*x11+x1*x3*x11+x2*x4*x11+x3*x5*x11+x4*x6*x11+x5*x7*x11+x6*x8*x11+x7*x9*x11+x8*x10*x11+x1*x11-2,
        x0*x3*x11+x1*x4*x11+x2*x5*x11+x3*x6*x11+x4*x7*x11+x5*x8*x11+x6*x9*x11+x7*x10*x11+x2*x11-3,
        x0*x4*x11+x1*x5*x11+x2*x6*x11+x3*x7*x11+x4*x8*x11+x5*x9*x11+x6*x10*x11+x3*x11-4,
        x0*x5*x11+x1*x6*x11+x2*x7*x11+x3*x8*x11+x4*x9*x11+x5*x10*x11+x4*x11-5,
        x0*x6*x11+x1*x7*x11+x2*x8*x11+x3*x9*x11+x4*x10*x11+x5*x11-6,
        x0*x7*x11+x1*x8*x11+x2*x9*x11+x3*x10*x11+x6*x11-7,
        x0*x8*x11+x1*x9*x11+x2*x10*x11+x7*x11-8,
        x0*x9*x11+x1*x10*x11+x8*x11-9,
        x0*x10*x11+x9*x11-10,
        x10*x11-11,
        x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+1
    ]
end

function eco13(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12) = np.PolynomialRing(ground, ["x$i" for i in 1:13], ordering=ordering)
    [
        x0*x1*x12+x1*x2*x12+x2*x3*x12+x3*x4*x12+x4*x5*x12+x5*x6*x12+x6*x7*x12+x7*x8*x12+x8*x9*x12+x9*x10*x12+x10*x11*x12+x0*x12-1,
        x0*x2*x12+x1*x3*x12+x2*x4*x12+x3*x5*x12+x4*x6*x12+x5*x7*x12+x6*x8*x12+x7*x9*x12+x8*x10*x12+x9*x11*x12+x1*x12-2,
        x0*x3*x12+x1*x4*x12+x2*x5*x12+x3*x6*x12+x4*x7*x12+x5*x8*x12+x6*x9*x12+x7*x10*x12+x8*x11*x12+x2*x12-3,
        x0*x4*x12+x1*x5*x12+x2*x6*x12+x3*x7*x12+x4*x8*x12+x5*x9*x12+x6*x10*x12+x7*x11*x12+x3*x12-4,
        x0*x5*x12+x1*x6*x12+x2*x7*x12+x3*x8*x12+x4*x9*x12+x5*x10*x12+x6*x11*x12+x4*x12-5,
        x0*x6*x12+x1*x7*x12+x2*x8*x12+x3*x9*x12+x4*x10*x12+x5*x11*x12+x5*x12-6,
        x0*x7*x12+x1*x8*x12+x2*x9*x12+x3*x10*x12+x4*x11*x12+x6*x12-7,
        x0*x8*x12+x1*x9*x12+x2*x10*x12+x3*x11*x12+x7*x12-8,
        x0*x9*x12+x1*x10*x12+x2*x11*x12+x8*x12-9,
        x0*x10*x12+x1*x11*x12+x9*x12-10,
        x0*x11*x12+x10*x12-11,
        x11*x12-12,
        x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+1
    ]
end

#------------------------------------------------------------------------------

function ku10(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = np.PolynomialRing(ground, ["x$i" for i in 1:10], ordering=ordering)
    [
        5*x1*x2+ 5*x1+ 3*x2+ 55,
        7*x2*x3+ 9*x2+ 9*x3+ 19,
        3*x3*x4+ 6*x3+ 5*x4-4,
        6*x4*x5+ 6*x4+ 7*x5+ 118,
        x5*x6+ 3*x5+ 9*x6+ 27,
        6*x6*x7+ 7*x6+x7+ 72,
        9*x7*x8+ 7*x7+x8+ 35,
        4*x8*x9+ 4*x8+ 6*x9+ 16,
        8*x9*x10+ 4*x9+ 3*x10-51,
        3*x1*x10-6*x1+x10+ 5
    ]
end

function kinema(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (z1, z2, z3, z4, z5, z6, z7, z8, z9) = np.PolynomialRing(ground, ["z$i" for i in 1:9], ordering=ordering)
    [
        z1^2 + z2^2 + z3^2 - 12*z1 - 68;
        z4^2 + z5^2 + z6^2 - 12*z5 - 68;
        z7^2 + z8^2 + z9^2 - 24*z8 - 12*z9 + 100;
        z1*z4 + z2*z5 + z3*z6 - 6*z1 - 6*z5 - 52;
        z1*z7 + z2*z8 + z3*z9 - 6*z1 - 12*z8 - 6*z9 + 64;
        z4*z7 + z5*z8 + z6*z9 - 6*z5 - 12*z8 - 6*z9 + 32;
        2*z2 + 2*z3 - z4 - z5 - 2*z6 - z7 - z9 + 18;
        z1 + z2 + 2*z3 + 2*z4 + 2*z6 - 2*z7 + z8 - z9 - 38;
        z1 + z3 - 2*z4 + z5 - z6 + 2*z7 - 2*z8 + 8;
    ]
end

function sparse5(; np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x1, x2, x3, x4, x5) = np.PolynomialRing(ground, ["x$i" for i in 1:5], ordering=ordering)
    [
        x1^2*x2^2*x3^2*x4^2*x5^2 + 3*x1^2 + x2^2 + x3^2 + x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + 3*x2^2 + x3^2 + x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x2^2 + 3*x3^2 + x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x2^2 + x3^2 + 3*x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x2^2 + x3^2 + x4^2 + 3*x5^2 + x1*x2*x3*x4*x5 + 5
    ]
end

function s9_1(; np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (a, b, c, d, e, f, g, h) = np.PolynomialRing(ground, ["x$i" for i in 1:8], ordering=ordering)
    [
        -e*g - 2*d*h,
        9*e + 4*b,
        -4*c*h - 2*e*f - 3*d*g,
        -7*c + 9*a - 8*f,
        -4*d*f - 5*c*g - 6*h - 3*e,
        -5*d - 6*c*f - 7*g + 9*b,
        9*d + 6*a - 5*b,
        9*c - 7*a + 8
    ]
end

function ojika4(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x1, x2, x3) = np.PolynomialRing(ground, ["x$i" for i in 1:3], ordering=ordering)
    [
        x1+x3*x1^3+x1*x3*x2^2-x1*x3,
        10*x2-2*x2*x3*x1^2-x3*x2^3-x2*x3,
        -6*x3^2*x1^4-3*x1^2*x2^2*x3^2-x3^2*x1^2+28*x3*x1^2 - 3*x3^2*x2^4+2*x3^2*x2^2+7*x3*x2^2+x3^2-11*x3+10
    ]
end

function ojika3_d1R2(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x1, x2, x3) = np.PolynomialRing(ground, ["x$i" for i in 1:3], ordering=ordering)
    [
        x1^3*x3+x1*x3*x2^2-x1*x3+x1,
        -2*x1^2*x3*x2-x3*x2^3-x3*x2+ 10*x2,
        -6*x1^4*x3^2-3*x1^2*x3^2*x2^2-3*x3^2*x2^4-x1^2*x3^2+ 2*x3^2*x2^2+ 28*x1^2*x3+ 7*x3*x2^2+x3^2-11*x3+ 10
    ]
end

function ojika4_d1R2_d2R5(;np=AbstractAlgebra, ground=np.QQ, ordering=:lex)
    _, (x1, x2, x3) = np.PolynomialRing(ground, ["x$i" for i in 1:3], ordering=ordering)
    [
        x1^3*x3+x1*x3*x2^2-x1*x3+x1,
        -2*x1^2*x3*x2-x3*x2^3-x3*x2+ 10*x2,
        -6*x1^4*x3^2-3*x1^2*x3^2*x2^2-3*x3^2*x2^4-x1^2*x3^2+ 2*x3^2*x2^2+ 28*x1^2*x3+ 7*x3*x2^2+x3^2-11*x3+ 10
    ]
end

#------------------------------------------------------------------------------

# generate a random polynomial system with:
# . `nvariables` variables, 
# . monomial exponents in range `exps`,
# . number of terms in range `nterms`,
# . `npolys` polynomials in total,
# . coefficients not greater than `csz`,
# . using `rng` random number generator,
# . over `ground` (must be a field of integers modulo prime),
# . in `ordering` term ordering
function generate_set(nvariables, exps, nterms, npolys, csz, rng,
                            ground, ordering; np=AbstractAlgebra)
    R, _ = np.PolynomialRing(
        ground,
        ["x$i" for i in 1:nvariables],
        ordering=ordering
    )

    filter!(!iszero, [
        np.map_coefficients(
            c -> ground(AbstractAlgebra.data(c) % csz),
            rand(rng, R, exps, nterms)
        )
        for _ in 1:rand(rng, npolys)
    ])
end

# -//-
# over `ground` (must be a field of rationals)
function generate_set(nvariables, exps, nterms, npolys, csz, rng,
                        ground::T, ordering; np=AbstractAlgebra) where {T<:AbstractAlgebra.Rationals}

    semiground = np.GF(2^31 - 1)
    R, _ = np.PolynomialRing(
        semiground,
        ["x$i" for i in 1:nvariables],
        ordering=ordering
    )

    csz = BigInt(csz)
    zzbase = BigInt(9223372036854775837)
    randzz() = rand(-zzbase:zzbase)

    filter!(!iszero, [
        np.map_coefficients(
            c -> ground(mod(randzz(), csz), mod(randzz(), csz) + 1),
            rand(rng, R, exps, nterms)
        )
        for _ in 1:rand(rng, npolys)
    ])
end

# -//-
# over `ground` (must be a ring of integers)
function generate_set(nvariables, exps, nterms, npolys, csz, rng,
        ground::T, ordering; np=AbstractAlgebra) where {T<:AbstractAlgebra.Integers}
    s = generate_set(nvariables, exps, nterms, npolys, csz, rng, QQ, ordering)
    s = map(f -> np.map_coefficients(c -> numerator(c), f), s)
    s = map(f -> np.change_coefficient_ring(np.ZZ, f), s)
    s
end

#------------------------------------------------------------------------------

change_ordering(f, ordering) = change_ordering([f], ordering)

# changes the ordering of set of polynomials `fs` into `ordering`
function change_ordering(fs::AbstractArray, ordering)
    R = parent(first(fs))
    Rord, _ = AbstractAlgebra.PolynomialRing(base_ring(R), string.(AbstractAlgebra.gens(R)), ordering=ordering)
    map(f -> AbstractAlgebra.change_base_ring(base_ring(R), f, parent=Rord), fs)
end
