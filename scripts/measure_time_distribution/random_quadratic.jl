using AbstractAlgebra, Groebner

function randpoly(x, d, coeffs=1:10)
    M = collect(monomials((sum(x)+1)^d))
    C = map(_ -> rand(coeffs), M)
    f = parent(x[1])(C, exponent_vector.(M, 1))
    f
end

function randsys(k, n, d, coeffs=1:10)
    R, x = polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=:degrevlex)
    [randpoly(x, d, coeffs) for _ in 1:n]
end

k = GF(2^30+3)
d = 2
for n in 2:12
    sys = randsys(k, n, d)
    println("n = $n:")
    @time gb = groebner(sys);
end

