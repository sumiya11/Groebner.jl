using AbstractAlgebra, Nemo
using Polymake

Nemo.height_bits(f::fmpq_mpoly) = maximum(Nemo.height_bits, collect(coefficients(f)))

newton_polytope(f) = newton_polytope([f])
function newton_polytope(I::Vector{T}) where {T}
    R = parent(first(I))
    I = deepcopy(I)
    prepend!(I, gens(R))
    prepend!(I, [one(R)])
    points = splat(vcat)(map(collect ∘ exponent_vectors, I))
    points = map(p -> prepend!(p, 1), points)
    points = reduce(vcat, map(transpose, points))
    polytope.Polytope(POINTS=points)
end

function effective_nullstellensatz_1(F)
    R = parent(first(F))
    s = length(F)
    n = nvars(R)
    d = maximum(Nemo.total_degree, F)
    h = maximum(Nemo.height_bits, F)
    n, d, h = map(BigInt, (n, d, h))
    @info "" n d h
    # bound on the degree
    D = 4 * n * d^n
    # bound on the height
    H = 4 * n * (n + 1) * d^n * (h + log(s) + (n + 7) * log(n + 1) * d)
    H = round(BigInt, H)
    (D=D, H=H)
end

function effective_nullstellensatz_2(F)
    R = parent(first(F))
    s = length(F)
    n = nvars(R)
    d = maximum(Nemo.total_degree, F)
    h = maximum(Nemo.height_bits, F)
    P = newton_polytope(F)
    V = Rational{BigInt}(P.VOLUME)
    n, d, h = map(BigInt, (n, d, h))
    @info "" n d h V
    # bound on the degree
    D = 2 * n^2 * d * V
    # bound on the height
    H = 2 * (n + 1)^3 * d * V * (h + log(s) + 2^(2n + 4) * d * log(d + 1))
    H = round(BigInt, H)
    (D=D, H=H)
end

R, (x1, x2) = Nemo.QQ["x1", "x2"]

I = [one(R), x1, x2, x1^2 * x2^2 - 1, x1 + x2]
I = Groebner.katsuran(2, np=Nemo) .^ 10
I = [x1]

P = newton_polytope(I)

effective_nullstellensatz_1(I)
effective_nullstellensatz_2(I)
