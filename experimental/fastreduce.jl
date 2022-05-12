









struct RewriteSystem
    tails  # vector of tails
    rules  # vector of rewriters
end

struct Rewriter
    idx::Integer
    monom
    terminal
end

function addrule!(RS, i, m)
    push!(RS.rules, Rewriter(i, m, isone(m)))
end

function addtail!(RS, t)
    push!(RS.tails, (leading_monomial(t), tail(t)))
    addrule!(RS, length(RS.tails), one(t))
end

function findrewriter(monom, RS)
    for rule in reverse(RS.rules)
        t, d = divides(monom, rule.monom*first(RS.tails[rule.idx]))
        if t
            @info "found rule" rule
            addrule!(RS, rule.idx, rule.monom*d)
            return d, rule.monom*last(RS.tails[rule.idx])
        end
    end
    return monom, zero(monom)
end

function rewrite(mult, f, RS)
    newmonoms = []
    @warn "rewriting" mult f
    for m in monomials(f)
        mm = mult*m
        mult1, rewriter = findrewriter(mm, RS)
        @info m mm mult1 rewriter
        if !iszero(rewriter)
            append!(newmonoms, .- rewrite(mult*mult1, rewriter, RS))
        else
            push!(newmonoms, mm)
        end
    end
    newmonoms
end

function rewrite(f, RS)
    sum(rewrite(one(f), f, RS))
end

"""


"""
function interreduce(gens::Vector{Poly}) where {Poly}
    sort!(gens, by=total_degree ∘ leading_monomial)

    RS = RewriteSystem([], [])

    # rewriters[1] = x1 + x2 + x3

    @warn "interreducing" gens

    newgens = similar(gens)
    newgens[1] = map_coefficients(c -> c // leading_coefficient(gens[1]), gens[1])

    addtail!(RS, newgens[1])

    i = 2

    while i <= length(gens)
        newgens[i] = rewrite(gens[i], RS)
        newgens[i] = map_coefficients(c -> c // leading_coefficient(newgens[i]), newgens[i])
        addtail!(RS, newgens[i])
        @error "End of $i" RS
        i += 1
    end

    newgens
end
