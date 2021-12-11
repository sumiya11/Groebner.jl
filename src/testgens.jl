
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

function eco5()
    R, (x1, x2, x3, x4, x5) = PolynomialRing(GF(1073741827), ["x$i" for i in 1:5])

    fs = [
    (x1 + x1*x2 + x2*x3 + x3*x4)*x5 - 1,
     (x2 + x1*x3 + x2*x4)*x5 - 2,
             (x3 + x1*x4)*x5 - 3,
                       x4*x5 - 4,
           x1 + x2 + x3 + x4 + 1
    ]
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
