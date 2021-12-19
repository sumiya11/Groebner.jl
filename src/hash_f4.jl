
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
mutable struct GroebnerState{Tv}
    G::Vector{HashedPolynomial{Tv}}
    P::Vector{Tuple{HashedPolynomial{Tv}, HashedPolynomial{Tv}}}
end

Base.iterate(x::GroebnerState) = (x.G, 1)
Base.iterate(x::GroebnerState, state) = state == 1 ? (x.P, 2) : nothing

GroebnerState{Tv}() where {Tv} = GroebnerState{Tv}([], [])

#------------------------------------------------------------------------------

function leading_monomials(polys::T) where {T<:AbstractArray}
    map(leading_monomial, polys)
end

function AbstractAlgebra.terms(polys::T) where {T<:AbstractArray}
    union(map(collect ∘ terms, polys)...)
end

# This is slow
function AbstractAlgebra.monomials(polys::Vector{T}) where {T}
    ans = Set{HashedMonomial}()
    for poly in polys
        for m in monomials(poly)
            push!(ans, m)
        end
    end
    collect(ans)
end

#------------------------------------------------------------------------------

# Returns degree of the given pair of polynomials, defined as
# the degree of the LCM of the leading terms of these polys
pairdegree(x, y, HT) = sum(HT[monom_lcm(leading_monomial(x), leading_monomial(y), HT)])

# Normal selection strategy
# Given an array of pairs, selects all of the lowest degree
# where the degree of (f, g) is equal to the degree of lcm(lm(f), lm(g))
function selectnormal(criticalpairs, HT; maxpairs=0)
    d = minimum(pairdegree(pair..., HT) for pair in criticalpairs)
    selected = filter(p -> pairdegree(p..., HT) == d, criticalpairs)

    if maxpairs != 0
        selected = selected[1:min(maxpairs, length(selected))]
    end

    filter!(p -> !(p in selected), criticalpairs)
    selected
end

#------------------------------------------------------------------------------

# Given pairs (f, g) constructs an array of pairs (t, f) such that
#   t*lead(f) == m*lead(g)
function leftright(
        ps::Vector{Tuple{HashedPolynomial{Tv}, HashedPolynomial{Tv}}},
        HT) where {Tv}

    left = Vector{Tuple{HashedMonomial, HashedPolynomial{Tv}}}(undef, 0)
    right = Vector{Tuple{HashedMonomial, HashedPolynomial{Tv}}}(undef, 0)

    for (fi, fj) in ps
        lti = leading_monomial(fi)
        ltj = leading_monomial(fj)
        lcmij = monom_lcm(lti, ltj, HT)

        push!(left, (monom_divide(lcmij, lti, HT), fi))
        push!(right, (monom_divide(lcmij, ltj, HT), fj))
    end

    union!(left, right)
end

#------------------------------------------------------------------------------

# constructs the matrix rows from the given polynomials F wrt monomials Tf
# and sorts them to obtain triangular shape
function constructrows_dense(F, Tf, HT)
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
function constructpolys_dense(rrefrows, Tf, HT)
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

# constructs the matrix rows from the given polynomials F wrt monomials Tf
# and sorts them to obtain triangular shape
#
# We returns the vectors sorted by pivot and then by sparsity
function constructrows_sparse(F, Tf, HT)
    # zero rows to fill
    ground = base_ring(first(F))
    rows = map(
        f -> SparseVectorAA(
                    length(Tf),
                    collect(1:length(F[f])),
                    zeros(ground, length(F[f])),
                    ground(0) ),
        1:length(F)
    )

    hash2idx = Dict(Tf[i] => i for i in 1:length(Tf))

    for i in 1:length(F)
        for (j, c, m) in zip(1:length(F[i]), coefficients(F[i]), monomials(F[i]))
            ind = hash2idx[m]
            # we can eliminate these allocations
            rows[i].nzval[j] = c
            rows[i].nzind[j] = ind
        end
    end

    triangular_pred = row -> (first(row.nzind), -nnz(row))
    sort!(rows, by=triangular_pred, rev=true)
end

# performs (almost backward) gaussian elimination
# in a dense deterministic format,
# assumes the given matrix is in triangular shape
function rref_sparse!(A)
    m, n = length(A), length(first(A))

    i = 1

    while i <= m
        while i <= m && nnz(A[i]) == 0
            i += 1
        end
        (i > m) && break

        pivot = first(nonzeroinds(A[i]))
        mul!(A[i], inv(A[i][pivot]))

        for j in i+1:m
            # iszero(A[j][pivot]) && break
            iszero(A[j][pivot]) && continue

            pivoto = A[j][pivot]
            mul!(A[i], pivoto)
            A[j] -= A[i]
            mul!(A[i], inv(pivoto))

            dropzeros!(A[j])
        end
        i += 1
    end

    A
end

# Constructs polynomials from the coefficients of rows
# given in rref
function constructpolys_sparse(rrefrows, Tf, targetring::MPolyRing{Tv}, HT) where {Tv}
    m, n = length(rrefrows), length(first(rrefrows))
    ground = base_ring(targetring)

    ans = Vector{HashedPolynomial{Tv}}(undef, 0)

    # for each row..
    for i in 1:m
        if nnz(rrefrows[i]) == 0
            continue
        end

        # we may eliminate sorting if the monomials are already ordered
        # important invariant: monomials are sorted in reverse order,
        #                     exponents in AA polys are stored in reverse order
        # another important invariant: all monomial orderings within this function
        #                             are same
        #
        poly = sparsevector_to_hashpoly(rrefrows[i], Tf, targetring, HT)
        push!(ans, poly)
    end

    ans
end

#------------------------------------------------------------------------------

# performs (almost backward) gaussian elimination
# in a dense deterministic format,
# assumes the given matrix is in triangular shape
function rref_sparse_randomized!(A)
    m, n = length(A), length(first(A))

    ground = base_ring(first(first(A)))

    ell = 5
    buf = zeros(ground, n)

    pivot_vectors = zeros(1)

    i = 1

    for sector in 1:m // ell
        for i in 1:ell
            for vector in sector
                cf = rand(ground)
                buf += cf * vector
            end
            reduce_vector!(buf, pivot_vectors)
            if !iszero(buf)
                add_vector!(pivot_vectors, buf)
            end
        end
        buf .= zero(ground)
    end

    A
end

#------------------------------------------------------------------------------

# constructs vectors from the given polynomials,
# then performs reduction (i.e, gaussian elimination),
# and transforms resulting row vectors back into polynomials
#
# This is standard dense linear algebra backend
function linear_algebra_dense_det(F, Tf, HT)
    A = constructrows_dense(F, Tf, HT)

    m, n = length(A), length(first(A))
    @debug "Constructed *dense* matrix of size $((m, n)) and $(sum(map(x->count(!iszero, x), A)) / (m*n)) nnz"

    Arref = rref_dense!(A)

    @debug "After reduction: $(sum(map(x->count(!iszero, x), Arref)) / (m*n)) nnz"

    constructpolys_dense(Arref, Tf, HT)
end

# Same as above but in sparse algebra
function linear_algebra_sparse_det(F, Tf, HT)
    A = constructrows_sparse(F, Tf, HT)

    m, n = length(A), length(first(A))
    @debug "Constructed *sparse* matrix of size $((m, n)) and $(sum(map(nnz, A)) / (m*n)) nnz"

    Arref = rref_sparse!(A)

    @debug "After reduction $(sum(map(nnz, A)) / (m*n)) nnz"

    constructpolys_sparse(Arref, Tf, parent(first(F)), HT)
end

# Same as above but in sparse random algebra
function linear_algebra_sparse_rand(F, Tf)
    A = constructrows_sparse(F, Tf)

    m, n = length(A), length(first(A))
    @debug "Constructed *sparse* matrix of size $((m, n)) and $(sum(map(nnz, A)) / (m*n)) nnz"

    Arref = rref_sparse_randomized!(A)

    constructpolys_sparse(Arref, Tf)
end

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

# Prepares the given set of polynomials L to be reduced by the ideal G,
# starts the reduction routine,
# and returns reduced polynomials
function reduction(L, G, HT; linalg=:sparse)

     @vtime F = symbolic_preprocessing(L, G, HT)

     cmp(x, y) = monom_is_less(x, y, parent(first(F)), HT)
     Tf = sort(monomials(F), lt=cmp, rev=true)

    @debug "Going to build a matrix of sizes $((length(F), length(Tf)))"

    if linalg == :dense
         F⁺ = linear_algebra_dense_det(F, Tf, HT)
    elseif linalg == :sparse
         @vtime F⁺ = linear_algebra_sparse_det(F, Tf, HT)
    elseif linalg == :randomsparse
         F⁺ = linear_algebra_sparse_rand(F, Tf)
    end

    Tf = Set(leading_monomials(F))

    F⁺ = filter!(f -> !in(leading_monomial(f), Tf), F⁺)

    F⁺, F
end

#------------------------------------------------------------------------------

# Given the set of polynomials L and the basis G,
# extends L to contain possible polys for reduction by G,
# and returns it
function symbolic_preprocessing(L, G, HT)
    PolyT = eltype(G)

    # F ~ several polys
    F     = unique(PolyT[monom_poly_product(t, f, HT) for (t, f) in L])
    Fhash = Set(F)

    Done = Set(leading_monomial(f) for f in F)

    # Tf ~ many monomials from F
    Tf = Set(filter(x -> !(x in Done), monomials(F))) # lm == buchberger

    # essentially for each in Tf
    while !isempty(Tf)
        m = pop!(Tf)
        # for each poly in G
        for f in G
            # PROFILER: this takes long
            # isdivisible: monom, lt(poly)
            flag = is_monom_divisible(m, leading_monomial(f), HT)
            if flag
                # divide: monom, lt(poly)
                t = term_divide(m, leading_term(f), HT)
                # multiplie: monom, poly
                tf = term_poly_product(t, f, HT)

                # poly in set of polys?
                if !(tf in Fhash)
                    push!(F,     tf)
                    push!(Fhash, tf)
                end

                # several monomials in set of monoms

                for (midx, ttf) in zip(1:length(tf), monomials(tf))
                    (midx == 1) && continue
                    push!(Tf, ttf)
                end
            end
        end
    end

    F
end

#------------------------------------------------------------------------------

# Hereinafter a set of heuristics is defined to be used in update! (see below)
# They assess which polynomials are worthy of adding to the current pairset

function update_heuristic_1(C, h, HT)
    D = similar(C, 0)
    while !isempty(C)
        h, g = pop!(C)

        flag1 = is_monom_gcd_constant(leading_monomial(g), leading_monomial(h), HT)
        flag2 = !any(t -> is_monom_divisible(
                            monom_lcm(leading_monomial(h), leading_monomial(g), HT),
                            monom_lcm(leading_monomial(h), leading_monomial(t[2]), HT),
                            HT),
                    C)
        flag3 = !any(t -> is_monom_divisible(
                            monom_lcm(leading_monomial(h), leading_monomial(g), HT),
                            monom_lcm(leading_monomial(h), leading_monomial(t[2]), HT),
                            HT),
                    D)

        if flag1 || (flag2 && flag3)
            push!(D, (h, g))
        end
    end
    D
end

function update_heuristic_2(D, h, HT)
    filter!(hg -> !is_monom_gcd_constant(leading_monomial(hg[1]), leading_monomial(hg[2]), HT), D)
end


function update_heuristic_3!(P, h, HT)
    flag1 = (g1, g2, h) -> !is_monom_divisible(monom_lcm(leading_monomial(g1), leading_monomial(g2), HT), leading_monomial(h), HT)
    flag2 = (g1, g2, h) -> monom_lcm(leading_monomial(g1), leading_monomial(h), HT) == monom_lcm(leading_monomial(g1), leading_monomial(g2), HT)
    flag3 = (g1, g2, h) -> monom_lcm(leading_monomial(g2), leading_monomial(h), HT) == monom_lcm(leading_monomial(g1), leading_monomial(g2), HT)
    flag  = (g1, g2, h) -> flag1(g1, g2, h) ||  flag2(g1, g2, h) || flag3(g1, g2, h)

    filter!(g -> flag(g[1], g[2], h), P)
end

function update_heuristic_4!(G, h, HT)
    filter!(g -> !is_monom_divisible(leading_monomial(g), leading_monomial(h), HT), G)
end

# "Adds" h to the set of generators G and set of pairs P, while applying some
# heuristics to reduce the number of pairs on fly
# Returns new generator and pair sets G, P
function update!(state, h, HT)
    # The algorithm and notation is taken from the book
    # Gröbner Bases - A Computational Approach to Commutative Algebra, 1993
    # Thomas Becker and Volker Weispfenning, Corrected second printing, page 230

    C = [(h, g) for g in state.G]

    # PROFILER: Usually the slowest one
    D = update_heuristic_1(C, h, HT)

    E = update_heuristic_2(D, h, HT)

    update_heuristic_3!(state.P, h, HT)
    append!(state.P, E)

    update_heuristic_4!(state.G, h, HT)
    push!(state.G, h)
end

#------------------------------------------------------------------------------

#
"""
    Main function to calculate the Groebner basis of the given polynomial ideal.
    Specialized to work only over finite fields.

    Parameters
        . F         - an array of polynomials over finite field
        . select    - a strategy for polynomial selection on each iteration
        . reduced   - reduce the basis so that the result is unique
        . linalg    - linear algebra backend to use. Possible options
                      are :dense , :sparse , and :sparserand. Default is :sparse
        . maxpairs  - maximal number of pairs selected for one matrix; default is
                      0, i.e. no restriction. If matrices get too big or consume
                      too much memory this is a good parameter to play with.
        . tablesize -
"""
function f4(polys::Vector{MPoly{GFElem{Int}}};
            select=selectnormal,
            reduced=true,
            linalg=:sparse,
            maxpairs=0,
            tablesize=2^16)

    # vector{MPoly}, vector{MPoly, MPoly}
    # stores currently accumulated groebner basis,
    # and an array of pairs to be reduced
    state = GroebnerState{GFElem{Int}}()

    HT = MonomialHashtable(parent(first(polys)))

    @vtime internal_polys = convert_to_hash_repr(polys, HT)

    # init with the given polys,
    # allow copy
    for poly in internal_polys
        update!(state, poly, HT)
    end

    d::Int64 = 0
    # while there are pairs to be reduced
    while !isempty(state.P)
        d += 1
        @debug "F4 iter $d"

        # vector{MPoly, MPoly}
        @vtime selected_pairs = select(state.P, HT, maxpairs=maxpairs)

        @debug "Seleted $(length(selected_pairs)) pairs for reduction"

        # vector{Monom, MPoly}
        @vtime to_be_reduced = leftright(selected_pairs, HT)

        # reduce polys and obtain new elements
        # to add to basis potentially
        @vtime Fd⁺, Fd = reduction(to_be_reduced, state.G, HT, linalg=linalg)

        # update the current basis,
        # allow copy
        @vtime for h in Fd⁺
            update!(state, h, HT)
        end

        if d > 1000
            @error "Something is probably wrong in f4.."
            break
        end
    end

    ans = convert_to_original_repr(state.G, HT)

    @debug "Filled $(length(HT.expmap)) elems in hashtable"

    if reduced
        @vtime ans = reducegb(ans)
    end

    ans
end

#------------------------------------------------------------------------------
