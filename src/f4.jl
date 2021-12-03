
#=
    The file contains implementation of the F4 algorithm,
    described originally by Jean-Charles Faugere in
        A new efficient algorithm for computing Grobner bases
=#

#=
    TODO:
        - warning if input set is empty / zero ?
=#

#------------------------------------------------------------------------------

# Stores the state of the Groebner basis and pairs to be assessed
struct GroebnerState
    G::Vector{MPoly{GFElem{Int}}}
    P::Vector{Tuple{MPoly{GFElem{Int}}, MPoly{GFElem{Int}}}}
end

Base.iterate(x::GroebnerState) = (x.G, 1)
Base.iterate(x::GroebnerState, state) = state == 1 ? (x.P, 2) : nothing

GroebnerState() = GroebnerState([], [])

#------------------------------------------------------------------------------

function leading_monomials(polys::T) where {T<:AbstractArray}
    map(leading_monomial, polys)
end

function AbstractAlgebra.terms(polys::T) where {T<:AbstractArray}
    union(map(collect ∘ terms, polys)...)
end

function AbstractAlgebra.monomials(polys::T) where {T<:AbstractArray}
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
# given polynomials F wrt monomials Tf
function constructmatrix(F, Tf)
    """
    Function is not used
    """
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
    """
    Function is not used
    """
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
function linear_algebra_aa(F)
    """
    Function is not used
    """
    Tf = sort(monomials(F), rev=true)

    MSpace = AbstractAlgebra.MatrixSpace(base_ring(first(F)), length(F), length(Tf))

    A = constructmatrix(F, Tf)

    Am = MSpace(A)
    Arref = rref(Am)

    constructpolys(Arref, Tf)
end

#------------------------------------------------------------------------------

# constructs the matrix rows from the given polynomials F wrt monomials Tf
# and sorts them to obtain triangular shape
function constructrows_dense(F, Tf)
    # zero rows to fill
    rows = map(_ -> zeros(base_ring(first(F)), length(Tf)), 1:length(F))
    for i in 1:length(F)
        for j in 1:length(Tf)
            rows[i][j] = coeff(F[i], Tf[j])
        end
    end
    triangular_pred = row -> (findfirst(!iszero, row), count(iszero, row))
    sort!(rows, by=triangular_pred, rev=true)
end

# performs (almost backward) gaussian elimination
# in a dense deterministic format,
# assumes the given matrix is in triangular shape
function rref_dense!(A)
    m, n = length(A), length(first(A))

    i = 1
    while i <= m
        while i <= m && iszero(A[i])
            i += 1
        end
        (i > m) && break

        pivot = findfirst(!iszero, A[i])
        A[i] .//= A[i][pivot]
        for j in i+1:m
            A[j] .-= A[i] .* A[j][pivot]
        end
        i += 1
    end

    A
end

# Constructs polynomials from the coefficients of rows
# given in rref
function constructpolys_dense(rrefrows, Tf)
    m, n = length(rrefrows), length(first(rrefrows))
    targetring = parent(first(Tf))
    ground = base_ring(targetring)

    ans = zeros(targetring, 0)

    # for each row..
    for i in 1:m
        if iszero(rrefrows[i])
            continue
        end

        builder = MPolyBuildCtx(targetring)
        # for each coeff in a row..
        for j in findfirst(!iszero, rrefrows[i]):n
            if iszero(rrefrows[i][j])
                continue
            end
            push_term!(builder, rrefrows[i][j], exponent_vector(Tf[j], 1))
        end

        push!(ans, finish(builder))
    end

    ans
end

#------------------------------------------------------------------------------

# constructs vectors from the given polynomials,
# then performs reduction (i.e, gaussian elimination),
# and transforms resulting row vectors back into polynomials
#
# This is standard dense linear algebra backend
function linear_algebra_dense_det(F)
    Tf = sort(monomials(F), rev=true)

    A = constructrows_dense(F, Tf)

    m, n = length(A), length(first(A))
    @debug "Constructed matrix of size $((m, n)) and $(sum(map(x->count(!iszero, x), A)) / (m*n)) nnz"

    Arref = rref_dense!(A)

    constructpolys_dense(Arref, Tf)
end

# Prepares the given set of polynomials L to be reduced by the ideal G,
# starts the reduction routine,
# and returns reduced polynomials
function reduction(L, G)
    F = symbolic_preprocessing(L, G)

    F⁺ = linear_algebra_dense_det(F)

    Tf = Set(leading_monomials(F))

    F⁺ = filter!(f -> !in(leading_monomial(f), Tf), F⁺)

    F⁺, F
end

#------------------------------------------------------------------------------

# Simplifies... well, not yet
function simplify(t, f)
    return t, f
end

# Given the set of polynomials L and the basis G,
# extends L to contain possible polys for reduction,
# and returns it
function symbolic_preprocessing(L, G)
    PolyT = eltype(G)

    F = unique(PolyT[t * f for (t, f) in L])
    Done = Set(leading_monomials(F))

    Tf = Set(filter(x -> !(x in Done), monomials(F)))

    while !isempty(Tf)
        m = pop!(Tf)
        push!(Done, m)
        for f in G
            flag, t = divides(m, leading_monomial(f))
            if flag
                tf = t * f
                Ttf = filter!(
                        x -> x != leading_monomial(tf),
                        collect(monomials(tf)))
                union!(F, [tf])
                union!(Tf, Ttf)
            end
        end
    end

    F
end

#------------------------------------------------------------------------------

# Hereinafter a set of heuristics is defined to be used in update! (see below)
# They assess which polynomials are worthy of adding to the current pairset

function update_heuristic_1(C, h)
    lmh = leading_monomial(h)
    D = similar(C, 0)
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
    D
end

function update_heuristic_2(D, h)
    lmh = leading_monomial(h)
    E = similar(D, 0)
    while !isempty(D)
        h, g = pop!(D)
        lmg = leading_monomial(g)
        if !isconstant(gcd(lmg, lmh))
            push!(E, (h, g))
        end
    end
    E
end

function update_heuristic_3!(E, h, P)
    lmh = leading_monomial(h)
    Pnew = similar(E, 0)
    while !isempty(P)
        g1, g2 = pop!(P)
        lmg1 = leading_monomial(g1)
        lmg2 = leading_monomial(g2)

        flag1 = !first(divides(lcm(lmg1, lmg2), lmh))
        flag2 = lcm(lmg1, lmh) == lcm(lmg1, lmg2)
        flag3 = lcm(lmg2, lmh) == lcm(lmg1, lmg2)
        if flag1 || flag2 || flag3
            push!(Pnew, (g1, g2))
        end
    end
    Pnew
end

function update_heuristic_4!(G, h)
    lmh = leading_monomial(h)
    Gnew = similar(G, 0)
    while !isempty(G)
        g = pop!(G)
        if !first(divides(leading_monomial(g), lmh))
            push!(Gnew, g)
        end
    end
    Gnew
end

# "Adds" h to the set of generators G and set of pairs P, while applying some
# heuristics to reduce the number of pairs on fly
# Returns new generator and pair sets G, P
function update!(state, h)
    # The algorithm and notation is taken from the book
    # Gröbner Bases - A Computational Approach to Commutative Algebra, 1993
    # Thomas Becker and Volker Weispfenning, Corrected second printing, page 230

    C = [(h, g) for g in state.G]

    D = update_heuristic_1(C, h)

    E = update_heuristic_2(D, h)

    Pnew = update_heuristic_3!(E, h, state.P)

    Gnew = update_heuristic_4!(state.G, h)

    union!(Pnew, E)
    push!(Gnew, h)

    copy!(state.G, Gnew)
    copy!(state.P, Pnew)
end

#------------------------------------------------------------------------------

#
"""
    Main function to calculate the Groebner basis of the given polynomial ideal.
    Specialized to work only over finite fields.

    Parameters
        . F       - an array of polynomials over finite field
        . select  - a strategy for polynomial selection on each iteration
        . reduced - reduce the basis so that the result is unique
"""
function f4(F::Vector{MPoly{GFElem{Int}}};
                select=selectnormal, reduced=true)

    # vector{MPoly}, vector{MPoly, MPoly}
    # stores currently accumulated groebner basis,
    # and an array of pairs to be reduced
    state = GroebnerState()

    # init with the given polys,
    # allow copy
    for f in F
        update!(state, f)
    end

    d = 0
    # while there are pairs to be reduced
    while !isempty(state.P)
        d += 1
        @debug "F4 iter $d"

        # vector{MPoly, MPoly}
        @vtime selected_pairs = select(state.P)

        # vector{Monom, MPoly}
        @vtime to_be_reduced = leftright(selected_pairs)

        # reduce polys and obtain new elements
        # to add to basis potentially
        @vtime Fd⁺, Fd = reduction(to_be_reduced, state.G)

        # update the current basis,
        # allow copy
        @vtime for h in Fd⁺
            update!(state, h)
        end

        if d > 100
            @error "Something is probably wrong in f4.."
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
