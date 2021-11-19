
"""
    The file contains implementation of the F4 algorithm,
    described originally by Jean-Charles Faugere in
        A new efficient algorithm for computing Grobner bases
"""

#=
    TODO:
        - warning if input set is empty / zero
=#

#------------------------------------------------------------------------------

# Stores the state of Groebner basis and pairs to assess
struct GroebnerState
    G::Vector{MPoly{GFElem{Int}}}
    P::Vector{Tuple{MPoly{GFElem{Int}}, MPoly{GFElem{Int}}}}
end

Base.iterate(x::GroebnerState) = (x.G, 1)
Base.iterate(x::GroebnerState, state) = state == 1 ? (x.P, 2) : nothing

GroebnerState() = GroebnerState([], [])

#------------------------------------------------------------------------------

function HT(polys)
    map(leading_monomial, polys)
end

function AbstractAlgebra.terms(polys)
    union(map(collect ∘ terms, polys)...)
end

function AbstractAlgebra.monomials(polys)
    union(map(collect ∘ monomials, polys)...)
end

#------------------------------------------------------------------------------

# Returns degree of the given pair of polynomials, defined as
# the degree of the LCM of the leading terms of these polys
pairdegree(x, y) = sum(degrees(lcm(leading_monomial(x), leading_monomial(y))))
pairdegree(p) = pairdegree(p...)

# Normal selection strategy
# Given an array of pairs, selects all of the lowest degree
# where the degree of (f, g) is equal to the degree of lcm(lm(f), lm(g))
function selectnormal(criticalpairs)
    d = minimum(pairdegree(pair) for pair in criticalpairs)
    selected = filter(p -> pairdegree(p) == d, criticalpairs)
    filter!(p -> !(p in selected), criticalpairs)
    selected
end

# Just fot testing:

# Strange selection strategy
# Given an array of pairs, selects and returns all
function selectall(polys)
    polys, []
end

# Strange selection strategy
# Given an array of pairs, selects and returns one
function selectone(polys)
    [polys[1]], polys[2:end]
end

#------------------------------------------------------------------------------

# Given pairs (f, g) constructs an array of pairs (t, f) such that
#   t*lead(f) == m*lead(g)
function leftright(ps)
    T = eltype(ps)

    left = T[]
    right = T[]

    for (fi, fj) in ps
        lti = leading_monomial(fi)
        ltj = leading_monomial(fj)
        lcmij = lcm(lti, ltj)

        push!(left, (div(lcmij, lti), fi))
        push!(right, (div(lcmij, ltj), fj))
    end

    union!(left, right)
end

#------------------------------------------------------------------------------

# Constructs a numeric matrix for reduction from the
# provided polynomials F and monomials Tf
function constructmatrix(F, Tf)
    shape = length(F), length(Tf)
    ground = base_ring(first(Tf))
    FF = elem_type(ground)

    A = Array{FF}(undef, shape...)
    for (i, f) in enumerate(F)
        for (j, t) in enumerate(Tf)
            A[i, j] = coeff(f, t)
        end
    end

    A
end

# Constructs polynomials from the rows coefficients of the given matrix Arref
# and the given monomials Tf
function constructpolys(Arref, Tf)
    k, A = Arref
    m, n = size(A)
    ground = base_ring(first(A))
    targetring = parent(first(Tf))

    ans = Array{elem_type(targetring), 1}()

    for i in 1:k
        f = targetring(0)
        for j in 1:n
            f += A[i, j] * Tf[j]
        end
        if f != 0
            push!(ans, f)
        end
    end

    ans
end

# Simultaneosly reduces the given set of polynomials by constructing a matrix
# and performing forward and backward gaussian elimination on it
function linear_algebra(F)
    Tf = sort(monomials(F), rev=true)

    MSpace = AbstractAlgebra.MatrixSpace(base_ring(first(F)), length(F), length(Tf))

    A = constructmatrix(F, Tf)
    @info "Matrix $(size(A))"
    if eltype(A) <: GFElem{Int}
        println(spy(AbstractAlgebra.lift.(A)))
    end
    Am = MSpace(A)
    Arref = rref(Am)

    constructpolys(Arref, Tf)
end

# Prepares the given set of polynomials L to be reduced by the ideal G
# and starts the reduction routine
function reduction(L, G, Fs)
    F = symbolic_preprocessing(L, G, Fs)

    F⁺ = linear_algebra(F)

    Tf = Set(HT(F))

    F⁺ = filter!(f -> !in(leading_monomial(f), Tf), F⁺)

    F⁺, F
end

#------------------------------------------------------------------------------

# Simplifies... well, not yet
function simplify(t, f, Fs)
    return t, f
end

# Eliminates redundant polynomials from L with respect to ideal G
function symbolic_preprocessing(L, G, Fs)
    PolyT = eltype(G)

    F = unique(PolyT[prod(simplify(t, f, Fs)) for (t, f) in L])
    Done = Set(HT(F))

    Tf = Set(filter(x -> !(x in Done), monomials(F)))

    while !isempty(Tf)
        m = pop!(Tf)
        push!(Done, m)
        for f in G
            flag, t = divides(m, leading_monomial(f))
            if flag
                tf = prod(simplify(t, f, Fs))
                Ttf = filter!(x -> x != leading_monomial(tf), collect(monomials(tf)))
                union!(F, [tf])
                union!(Tf, Ttf)
            end
        end
    end

    F
end

#------------------------------------------------------------------------------

# "Adds" h to the set of generators G and set of pairs P, while applying some
# heuristics to reduce the number of pairs on fly
# Returns new generator and pair sets G, P
function update!(state, h)
    # The algorithm and notation is taken from the book
    # Gröbner Bases - A Computational Approach to Commutative Algebra, 1993
    # Thomas Becker and Volker Weispfenning, Corrected second printing, page 230

    G, P = state
    PolyT = eltype(G)

    C = Tuple{PolyT, PolyT}[(h, g) for g in G]
    D = similar(C, 0)

    lmh = leading_monomial(h)

    # heuristic 1
    while !isempty(C)
        h, g = pop!(C)
        lmg = leading_monomial(g)

        flag1 = isconstant(gcd(lmg, lmh))
        flag2 = all(t -> !first(divides(
                            lcm(lmh, lmg),
                            lcm(lmh, leading_monomial(t[2])))),
                    C)
        flag3 = all(t -> !first(divides(
                            lcm(lmh, lmg),
                            lcm(lmh, leading_monomial(t[2])))),
                    D)

        if flag1 || (flag2 && flag3)
            push!(D, (h, g))
        end
    end

    # heuristic 2
    E = similar(D, 0)
    while !isempty(D)
        h, g = pop!(D)
        lmg = leading_monomial(g)
        if !isconstant(gcd(lmg, lmh))
            push!(E, (h, g))
        end
    end

    # heuristic 3
    Pnew = similar(E, 0)
    while !isempty(state.P)
        g1, g2 = pop!(state.P)
        lmg1 = leading_monomial(g1)
        lmg2 = leading_monomial(g2)

        flag1 = !first(divides(lcm(lmg1, lmg2), lmh))
        flag2 = lcm(lmg1, lmh) == lcm(lmg1, lmg2)
        flag3 = lcm(lmg2, lmh) == lcm(lmg1, lmg2)
        if flag1 || flag2 || flag3
            push!(Pnew, (g1, g2))
        end
    end

    union!(Pnew, E)

    # heuristic 4
    Gnew = similar(state.G, 0)
    while !isempty(state.G)
        g = pop!(state.G)
        if !first(divides(leading_monomial(g), lmh))
            push!(Gnew, g)
        end
    end

    push!(Gnew, h)

    copy!(state.G, Gnew)
    copy!(state.P, Pnew)
end

#------------------------------------------------------------------------------

#
# Main function to calculate the Groebner basis of the given polynomial ideal
#
# F      : an array of polynomials
# select : strategy for selection to perform on each iteration.
#          Must accept an array of pairs and return its partition (selected, not)
# reduced: reduce the basis so that the result is unique
function f4(F::Vector{MPoly{GFElem{Int}}};
            select=selectnormal, reduced=true)

    state = GroebnerState()
    for f in F
        update!(state, f)
    end

    Fs = []

    d = 0
    while !isempty(state.P)
        Pd = select(state.P)

        Ld = leftright(Pd)

        Fd⁺, Fd = reduction(Ld, state.G, Fs)
        # push!(Fs, Fd)

        for h in Fd⁺
            update!(state, h)
        end

        d += 1
        @info "d = $d"

        if d > 100
            @error "Something is probably wrong.."
            break
        end
    end

    if reduced
        reduced_gb = reducegb(state.G)
        copy!(state.G, reduced_gb)
    end

    state.G
end

#------------------------------------------------------------------------------

function groebner(
    fs::Vector{MPoly{Rational{BigInt}}})

    

end

function groebner(fs::Vector{MPoly{GFElem{Int}}})
    f4(fs)
end

#------------------------------------------------------------------------------

#=
ground = GF(2^31 - 1)

gb = f4(rootn(5, ground=ground))

println(gb)


R, (x, y, z, w) = PolynomialRing(ground, ["x", "y", "z", "w"])
fs = [
    x + y,
    x^2 + z^2,
    y^2 + w^2,
    x + 2*y + 3*z + 4*w
]
G = f4(fs, reduced=true)

println(G)
=#
