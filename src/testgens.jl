
# nearest TODO: generate answers for these

function rootn(n; ground=QQ)
    R, xs = PolynomialRing(ground, ["x$i" for i in 1:n])
    ans = [
        sum(map(prod, combinations(xs, i)))
        for i in 1:n
    ]
    ans[end] -= (-1)^(n - 1)
    ans
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

function noon4()
    R, (x1, x2, x3, x4) = PolynomialRing(GF(1073741827), ["x$i" for i in 1:4])

    fs = [
    10*x1^2*x4+10*x2^2*x4+10*x3^2*x4-11*x4+10,
    10*x1^2*x3+10*x2^2*x3+10*x3*x4^2-11*x3+10,
    10*x1*x2^2+10*x1*x3^2+10*x1*x4^2-11*x1+10,
    10*x1^2*x2+10*x2*x3^2+10*x2*x4^2-11*x2+10
    ]
    fs
end

function noon5()
    R, (x1, x2, x3, x4, x5) = PolynomialRing(GF(1073741827), ["x$i" for i in 1:5])

    fs = [
    10*x1^2*x5+10*x2^2*x5+10*x3^2*x5+10*x4^2*x5-11*x5+10,
    10*x1^2*x4+10*x2^2*x4+10*x3^2*x4+10*x4*x5^2-11*x4+10,
    10*x1^2*x3+10*x2^2*x3+10*x3*x4^2+10*x3*x5^2-11*x3+10,
    10*x1*x2^2+10*x1*x3^2+10*x1*x4^2+10*x1*x5^2-11*x1+10,
    10*x1^2*x2+10*x2*x3^2+10*x2*x4^2+10*x2*x5^2-11*x2+10
    ]
    fs
end

function noon6()
    R, (x1, x2, x3, x4, x5, x6) = PolynomialRing(GF(1073741827), ["x$i" for i in 1:6])

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
