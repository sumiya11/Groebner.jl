
using Nemo
using Groebner
using Logging

Logging.global_logger(ConsoleLogger(Logging.Info))

function check_bases_shape(groebner_bases)
    @assert !isempty(groebner_bases)
    ell = length(first(groebner_bases))
    ka = map(length, first(groebner_bases))
    @assert all(gb -> length(gb) == ell, groebner_bases)
    @assert all(gb -> map(length, gb) == ka, groebner_bases)
end

# Groebner bases over the field of rational functions
# over a field
function groebner_formal_parametric(F)
    R = parent(first(F))
    R1 = base_ring(first(F))
    ps = gens(R1)
    np = length(ps)

    # univariatize the ring of parameters
    R1u, u = polynomial_ring(base_ring(R1), "u")
    # maximal degree
    d = 10
    seq = [d^(i - 1) for i in 1:np]
    us = [u^seq[i] for i in 1:np]
    Fu = map(f -> map_coefficients(c -> evaluate(c, us), f), F)
    @info "univariate" Fu

    # Generating generic points
    generic_points = [rand(1:(2^20)) for _ in 1:50]
    @info "points" generic_points

    # Computing specialized groebner bases
    groebner_bases = []
    for point in generic_points
        Fp = map(f -> map_coefficients(c -> evaluate(c, point), f), Fu)
        push!(groebner_bases, groebner(Fp))
    end
    @info "bases" groebner_bases

    # check bases correctness
    check_bases_shape(groebner_bases)

    # interpolate coefficients
    basis = [zeros(R1u, length(f)) for f in first(groebner_bases)]
    for (i, f) in enumerate(basis)
        for (j, c) in enumerate(f)
            xs = [QQ(point[1]) for point in generic_points]
            ys = [coeff(gb[i], j) for gb in groebner_bases]

            @info "interpol" i j
            println(xs, " --> ", ys)

            c = interpolate(R1u, xs, ys)

            @info "=" c

            basis[i][j] = c
        end
    end
    @info "coeffs" basis

    # multivariaze the ring of parameters back
    invseq = map(f -> map(c -> map(cc -> digits(cc, base=d), 0:degree(c)), f), basis)
    basis = [
        [
            sum(
                map(
                    prod,
                    zip(coefficients(basis[i][j]), map(ss -> prod(ps .^ ss), invseq[i][j]))
                )
            ) for j in 1:length(basis[i])
        ] for i in 1:length(basis)
    ]
    @info "multivariate" basis

    # construct the answer
    groebner_bases =
        map(gb -> map(f -> change_base_ring(R1, f, parent=R), gb), groebner_bases)
    [
        sum(map(prod, zip(basis[i], monomials(f)))) for
        (i, f) in enumerate(first(groebner_bases))
    ]
end

R1, (a, b) = polynomial_ring(QQ, ["a", "b"])
R, (x, y) = R1["x", "y"]

f1 = x + a + b
f2 = a * b * y

F = [f1, f2]

groebner_formal_parametric(F)
