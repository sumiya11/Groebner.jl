
# The file contains test examples definitions

# nearest TODO: generate answers for these
# and also consider doing this in another way

function rootn(n; ground=QQ)
    R, xs = PolynomialRing(ground, ["x$i" for i in 1:n])
    ans = [
        sum(map(prod, Combinatorics.combinations(xs, i)))
        for i in 1:n
    ]
    ans[end] -= (-1)^(n - 1)
    ans
end

function henrion5(;ground=QQ)

    R, (f1,f2,f3,f4,f5,t) = PolynomialRing(ground, ["f1","f2","f3","f4","f5","t"])
    fs = [
        2*f1*f2*f3*f4*f5-9823275,
        21//5*f1*f2*f4*f5+16//5*f1*f3*f4*f5+9//5*f2*f3*f4*f5+24//5*f1*f2*f3*f5+5*f4*f3*f1*f2-4465125,
        14//5*f4*f5*f1+14//5*f4*f5*f2+8//5*f3*f4*f5+18//5*f1*f2*f5+24//5*f1*f3*f5+18//5*f2*f3*f5+4*f3*f1*f2+6*f1*f2*f4+6*f3*f4*f1+4*f2*f3*f4-441486,
        7//5*f4*f5+12//5*f5*f1+12//5*f5*f2+12//5*f5*f3+3*f1*f2+4*f3*f1+4*f4*f1+3*f2*f3+4*f4*f2+3*f3*f4-15498,
        6//5*f5+2*f4+2*f3+2*f2+2*f1-215,
        f1+2*f2+3*f3+4*f4+5*f5+6*t
    ]
end

function katsura6(;ground=QQ)
    R, (x1, x2, x3, x4, x5, x6, x7) = PolynomialRing(ground, ["x$i" for i in 1:7])

    fs = [
        1*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7-1,
        2*x4*x3+2*x5*x2+2*x6*x1+2*x7*x2-1*x6,
        1*x3^2+2*x4*x2+2*x5*x1+2*x6*x2+2*x7*x3-1*x5,
        2*x3*x2+2*x4*x1+2*x5*x2+2*x6*x3+2*x7*x4-1*x4,
        1*x2^2+2*x3*x1+2*x4*x2+2*x5*x3+2*x6*x4+2*x7*x5-1*x3,
        2*x2*x1+2*x3*x2+2*x4*x3+2*x5*x4+2*x6*x5+2*x7*x6-1*x2,
        1*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2-1*x1
    ]
    fs
end

function katsura9(;ground=QQ)
    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9) = PolynomialRing(ground, ["x$i" for i in 1:9])

    fs = [
        x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9-1,
        x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2-x1,
        2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9-x2,
        x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9-x3,
        2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9-x4,
        x3^2+2*x2*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9-x5,
        2*x3*x4+2*x2*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9-x6,
        x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x2*x8+2*x3*x9-x7,
        2*x4*x5+2*x3*x6+2*x2*x7+2*x1*x8+2*x2*x9-x8
    ]
end

function katsura10(;ground=QQ)
    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = PolynomialRing(ground, ["x$i" for i in 1:10])

    fs = [
        x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9+2*x10-1,
        x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2+2*x10^2-x1,
        2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9+2*x9*x10-x2,
        x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9+2*x8*x10-x3,
        2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9+2*x7*x10-x4,
        x3^2+2*x2*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9+2*x6*x10-x5,
        2*x3*x4+2*x2*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9+2*x5*x10-x6,
        x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x2*x8+2*x3*x9+2*x4*x10-x7,
        2*x4*x5+2*x3*x6+2*x2*x7+2*x1*x8+2*x2*x9+2*x3*x10-x8,
        x5^2+2*x4*x6+2*x3*x7+2*x2*x8+2*x1*x9+2*x2*x10-x9
    ]
end

function katsura11(;ground=QQ)
    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) = PolynomialRing(ground, ["x$i" for i in 1:11])

    fs = [
        x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9+2*x10+2*x11-1,
        x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2+2*x10^2+2*x11^2-x1,
        2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9+2*x9*x10+2*x10*x11-x2,
        x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9+2*x8*x10+2*x9*x11-x3,
        2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9+2*x7*x10+2*x8*x11-x4,
        x3^2+2*x2*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9+2*x6*x10+2*x7*x11-x5,
        2*x3*x4+2*x2*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9+2*x5*x10+2*x6*x11-x6,
        x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x2*x8+2*x3*x9+2*x4*x10+2*x5*x11-x7,
        2*x4*x5+2*x3*x6+2*x2*x7+2*x1*x8+2*x2*x9+2*x3*x10+2*x4*x11-x8,
        x5^2+2*x4*x6+2*x3*x7+2*x2*x8+2*x1*x9+2*x2*x10+2*x3*x11-x9,
        2*x5*x6+2*x4*x7+2*x3*x8+2*x2*x9+2*x1*x10+2*x2*x11-x10
    ]
end

function katsura12(;ground=QQ)
    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12) = PolynomialRing(ground, ["x$i" for i in 1:12])

    fs = [
        x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7+2*x8+2*x9+2*x10+2*x11+2*x12-1,
        x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2+2*x10^2+2*x11^2+2*x12^2-x1,
        2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9+2*x9*x10+2*x10*x11+2*x11*x12-x2,
        x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9+2*x8*x10+2*x9*x11+2*x10*x12-x3,
        2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9+2*x7*x10+2*x8*x11+2*x9*x12-x4,
        x3^2+2*x2*x4+2*x1*x5+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9+2*x6*x10+2*x7*x11+2*x8*x12-x5,
        2*x3*x4+2*x2*x5+2*x1*x6+2*x2*x7+2*x3*x8+2*x4*x9+2*x5*x10+2*x6*x11+2*x7*x12-x6,
        x4^2+2*x3*x5+2*x2*x6+2*x1*x7+2*x2*x8+2*x3*x9+2*x4*x10+2*x5*x11+2*x6*x12-x7,
        2*x4*x5+2*x3*x6+2*x2*x7+2*x1*x8+2*x2*x9+2*x3*x10+2*x4*x11+2*x5*x12-x8,
        x5^2+2*x4*x6+2*x3*x7+2*x2*x8+2*x1*x9+2*x2*x10+2*x3*x11+2*x4*x12-x9,
        2*x5*x6+2*x4*x7+2*x3*x8+2*x2*x9+2*x1*x10+2*x2*x11+2*x3*x12-x10,
        x6^2+2*x5*x7+2*x4*x8+2*x3*x9+2*x2*x10+2*x1*x11+2*x2*x12-x11
    ]
end

function eco5(;ground=GF(2^31-1))
    R, (x1, x2, x3, x4, x5) = PolynomialRing(ground, ["x$i" for i in 1:5])

    fs = [
    (x1 + x1*x2 + x2*x3 + x3*x4)*x5 - 1,
     (x2 + x1*x3 + x2*x4)*x5 - 2,
             (x3 + x1*x4)*x5 - 3,
                       x4*x5 - 4,
           x1 + x2 + x3 + x4 + 1
    ]
end

function eco7(;ground=GF(2^31-1))
    R, (x1, x2, x3, x4, x5, x6, x7) = PolynomialRing(ground, ["x$i" for i in 1:7])

    fs = [
        (x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6)*x7 - 1,
     (x2 + x1*x3 + x2*x4 + x3*x5 + x4*x6)*x7 - 2,
     (x3 + x1*x4 + x2*x5 + x3*x6)*x7 - 3,
     (x4 + x1*x5 + x2*x6)*x7 - 4,
     (x5 + x1*x6)*x7 - 5,
     x6*x7 - 6,
     x1 + x2 + x3 + x4 + x5 + x6 + 1
    ]
    fs
end

function eco10(;ground=GF(2^31-1))
    R, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9) = PolynomialRing(ground, ["x$i" for i in 1:11])

    fs = [
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
    fs
end

function eco11(;ground=GF(2^31-1))
    R, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = PolynomialRing(ground, ["x$i" for i in 1:11])

    fs = [
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
    fs
end

function eco12(;ground=GF(2^31-1))
    R, (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) = PolynomialRing(ground, ["x$i" for i in 1:11])

    fs = [
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
    fs
end

function noon3(;ground=QQ)
    R, (x1, x2, x3) = PolynomialRing(ground, ["x$i" for i in 1:3])
    fs = [
    10x1*x2^2 + 10x1*x3^2 - 11x1 + 10,
    10x2*x1^2 + 10x2*x3^2 - 11x2 + 10,
    10x3*x1^2 + 10x3*x2^2 - 11x3 + 10,
    ]
    fs
end

function noon4(; ground=QQ)
    R, (x1, x2, x3, x4) = PolynomialRing(ground, ["x$i" for i in 1:4])

    fs = [
    10*x1^2*x4+10*x2^2*x4+10*x3^2*x4-11*x4+10,
    10*x1^2*x3+10*x2^2*x3+10*x3*x4^2-11*x3+10,
    10*x1*x2^2+10*x1*x3^2+10*x1*x4^2-11*x1+10,
    10*x1^2*x2+10*x2*x3^2+10*x2*x4^2-11*x2+10
    ]
    fs
end

function noon5(;ground=QQ)
    R, (x1, x2, x3, x4, x5) = PolynomialRing(ground, ["x$i" for i in 1:5])

    fs = [
    10*x1^2*x5+10*x2^2*x5+10*x3^2*x5+10*x4^2*x5-11*x5+10,
    10*x1^2*x4+10*x2^2*x4+10*x3^2*x4+10*x4*x5^2-11*x4+10,
    10*x1^2*x3+10*x2^2*x3+10*x3*x4^2+10*x3*x5^2-11*x3+10,
    10*x1*x2^2+10*x1*x3^2+10*x1*x4^2+10*x1*x5^2-11*x1+10,
    10*x1^2*x2+10*x2*x3^2+10*x2*x4^2+10*x2*x5^2-11*x2+10
    ]
    fs
end

function noon6(;ground=QQ)
    R, (x1, x2, x3, x4, x5, x6) = PolynomialRing(ground, ["x$i" for i in 1:6])

    fs = [
    10*x1^2*x6+10*x2^2*x6+10*x3^2*x6+10*x4^2*x6+10*x5^2*x6-11*x6+10,
    10*x1^2*x5+10*x2^2*x5+10*x3^2*x5+10*x4^2*x5+10*x5*x6^2-11*x5+10,
    10*x1^2*x4+10*x2^2*x4+10*x3^2*x4+10*x4*x5^2+10*x4*x6^2-11*x4+10,
    10*x1^2*x3+10*x2^2*x3+10*x3*x4^2+10*x3*x5^2+10*x3*x6^2-11*x3+10,
    10*x1*x2^2+10*x1*x3^2+10*x1*x4^2+10*x1*x5^2+10*x1*x6^2-11*x1+10,
    10*x1^2*x2+10*x2*x3^2+10*x2*x4^2+10*x2*x5^2+10*x2*x6^2-11*x2+10
    ]
    fs
end

function noonn(n; ground=QQ)
    without(x, k) = x[1:end .!= k]

    R, xs = PolynomialRing(ground, ["x$i" for i in 1:n])
    fs = zeros(R, n)
    for i in 1:n
        other = without(xs, i)
        fs[i] = xs[i] * (10*sum(other .^ 2) - 11) + 10
    end
    fs
end

function ku10(;ground=QQ)
    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) = PolynomialRing(ground, ["x$i" for i in 1:10])

    fs = [
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
    fs
end

function kinema(;ground=QQ)
    R, (z1, z2, z3, z4, z5, z6, z7, z8, z9) = PolynomialRing(ground, ["z$i" for i in 1:9])

    fs = [
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

function sparse5(; ground=QQ)
    R, (x1, x2, x3, x4, x5) = PolynomialRing(ground, ["x$i" for i in 1:5])

    fs = [
        x1^2*x2^2*x3^2*x4^2*x5^2 + 3*x1^2 + x2^2 + x3^2 + x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + 3*x2^2 + x3^2 + x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x2^2 + 3*x3^2 + x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x2^2 + x3^2 + 3*x4^2 + x5^2 + x1*x2*x3*x4*x5 + 5,
        x1^2*x2^2*x3^2*x4^2*x5^2 + x1^2 + x2^2 + x3^2 + x4^2 + 3*x5^2 + x1*x2*x3*x4*x5 + 5
    ]
end

function s9_1(; ground=QQ)
    R, (a, b, c, d, e, f, g, h) = PolynomialRing(ground, ["x$i" for i in 1:8])

    fs = [-e*g - 2*d*h,
        9*e + 4*b,
        -4*c*h - 2*e*f - 3*d*g,
        -7*c + 9*a - 8*f,
        -4*d*f - 5*c*g - 6*h - 3*e,
        -5*d - 6*c*f - 7*g + 9*b,
        9*d + 6*a - 5*b,
        9*c - 7*a + 8]
end

function ojika4(;ground=QQ)
    R, (x1, x2, x3) = PolynomialRing(ground, ["x$i" for i in 1:3])

    fs = [x1+x3*x1^3+x1*x3*x2^2-x1*x3,
    10*x2-2*x2*x3*x1^2-x3*x2^3-x2*x3,
    -6*x3^2*x1^4-3*x1^2*x2^2*x3^2-x3^2*x1^2+28*x3*x1^2 - 3*x3^2*x2^4+2*x3^2*x2^2+7*x3*x2^2+x3^2-11*x3+10
    ]
end

function ojika3_d1R2(;ground=QQ)
    R, (x1, x2, x3) = PolynomialRing(ground, ["x$i" for i in 1:3])

    fs = [x1^3*x3+x1*x3*x2^2-x1*x3+x1,
-2*x1^2*x3*x2-x3*x2^3-x3*x2+ 10*x2,
-6*x1^4*x3^2-3*x1^2*x3^2*x2^2-3*x3^2*x2^4-x1^2*x3^2+ 2*x3^2*x2^2+ 28*x1^2*x3+ 7*x3*x2^2+x3^2-11*x3+ 10]
end

function ojika4_d1R2_d2R5(;ground=QQ)
    R, (x1, x2, x3) = PolynomialRing(ground, ["x$i" for i in 1:3])

    fs = [
    x1^3*x3+x1*x3*x2^2-x1*x3+x1,
-2*x1^2*x3*x2-x3*x2^3-x3*x2+ 10*x2,
-6*x1^4*x3^2-3*x1^2*x3^2*x2^2-3*x3^2*x2^4-x1^2*x3^2+ 2*x3^2*x2^2+ 28*x1^2*x3+ 7*x3*x2^2+x3^2-11*x3+ 10
    ]
end

function generate_set(nvariables, exps, nterms, npolys, csz, rng,
                            ground, ordering)
    R, _ = PolynomialRing(
        ground,
        ["x$i" for i in 1:nvariables],
        ordering=ordering
    )

    return filter!(!iszero, [
        map_coefficients(
            c -> ground(data(c) % csz),
            rand(rng, R, exps, nterms)
        )
        for _ in 1:rand(rng, npolys)
    ])
end

function generate_set(nvariables, exps, nterms, npolys, csz, rng,
                        ground::AbstractAlgebra.Rationals, ordering)

    semiground = GF(2^31 - 1)
    R, _ = PolynomialRing(
        semiground,
        ["x$i" for i in 1:nvariables],
        ordering=ordering
    )

    csz = BigInt(csz)
    zzbase = BigInt(9223372036854775837)
    randzz() = rand(-zzbase:zzbase)

    return filter!(!iszero, [
        map_coefficients(
            c -> ground(mod(randzz(), csz), mod(randzz(), csz) + 1),
            rand(rng, R, exps, nterms)
        )
        for _ in 1:rand(rng, npolys)
    ])
end
