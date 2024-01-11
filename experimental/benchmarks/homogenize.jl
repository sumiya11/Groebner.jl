function homogenize(fs)
    ring = parent(fs[1])
    newring, hom_vars = polynomial_ring(
        base_ring(ring),
        vcat("X0", map(string, gens(ring))),
        ordering=ordering(ring)
    )
    Fs = empty(fs)
    for f in fs
        D = total_degree(f)
        new_f = zero(newring)
        for term in terms(f)
            cf = coeff(term, 1)
            ev = monomial(term, 1)
            d = total_degree(ev)
            new_f += cf * evaluate(ev, hom_vars[2:end]) * hom_vars[1]^(D - d)
        end
        push!(Fs, new_f)
    end
    return Fs
end

function dehomogenize(Fs)
    ring = parent(Fs[1])
    newring, dehom_vars = polynomial_ring(
        base_ring(ring),
        map(string, gens(ring)[2:end]),
        ordering=ordering(ring)
    )
    fs = empty(Fs)
    for F in Fs
        f = evaluate(F, vcat(one(newring), dehom_vars))
        push!(fs, f)
    end
    return fs
end

begin
    using AbstractAlgebra
    K = GF(2^31 - 1)
    kat = Groebner.katsuran(6, k=K)

    kat_hom = homogenize(kat)
    xs = gens(parent(kat_hom[1]))
    new_ring, new_xs = polynomial_ring(K, vcat("t", map(string, xs)))

    kat_hom_true = map(f -> evaluate(f, new_xs[2:end]), kat_hom)
    kat_hom_true = vcat(kat_hom_true, new_xs[1] * new_xs[2] - 1)
end

@time Groebner.groebner(kat, ordering=Groebner.DegRevLex());

@time Groebner.groebner(kat_hom, ordering=Groebner.DegRevLex());

@time Groebner.groebner(kat_hom_true, ordering=Groebner.DegRevLex());

@time Groebner.groebner(kat_hom, ordering=Groebner.Lex(), loglevel=-3);
@time Groebner.groebner(kat_hom_true, ordering=Groebner.Lex());
