# FGLM implementation is correct, but is still WIP

mutable struct NextMonomials{Ord}
    # monomials to check
    monoms::Vector{MonomIdx}
    load::Int
    done::Dict{Int, Int}
    ord::Ord
end

Base.isempty(m::NextMonomials) = m.load == 0

function initialize_nextmonomials(ht::MonomialHashtable{M}, ord) where {M}
    zz = construct_const_monom(M, ht.nvars)
    vidx = insert_in_hash_table!(ht, zz)
    monoms = Vector{MonomIdx}(undef, 2^3)
    monoms[1] = vidx
    load = 1
    NextMonomials{typeof(ord)}(monoms, load, Dict{MonomIdx, Int}(vidx => 1), ord)
end

function insertnexts!(m::NextMonomials, ht::MonomialHashtable{M}, monom::MonomIdx) where {M}
    while m.load + ht.nvars >= length(m.monoms)
        resize!(m.monoms, length(m.monoms) * 2)
    end

    emonom = ht.monoms[monom]
    for i in 1:(ht.nvars)
        eprod = copy(emonom)
        tmp = Vector{UInt}(undef, ht.nvars)
        edense = monom_to_dense_vector!(tmp, eprod)
        edense[i] += 1
        eprod = construct_monom(M, edense)
        vidx = insert_in_hash_table!(ht, eprod)

        if !haskey(m.done, vidx)
            m.load += 1
            m.monoms[m.load] = vidx
            m.done[vidx] = 1
        end
    end

    sort_monom_indices_decreasing!(m.monoms, m.load, ht, m.ord)
end

function nextmonomial!(m::NextMonomials)
    # assuming m.monoms is sorted (reversed)
    monom = m.monoms[m.load]
    m.load -= 1
    monom
end

function add_generator!(basis::Basis{C}, matrix, relation, ht, ord) where {C <: Coeff}
    rexps, rcoeffs, _ = extract_sparse_row(relation)

    # if debug()
    #     # @warn "add generator"
    #     println(relation)
    #     println(rexps, " ", rcoeffs)
    # end

    for i in 1:length(rexps)
        rexps[i] = matrix.rightcol2hash[rexps[i]]
    end

    sort_term_indices_decreasing!(rexps, rcoeffs, ht, ord)

    # if debug()
    #     # @warn "extracted"
    #     println(rexps)
    #     println(rcoeffs)
    # end

    resize_basis_if_needed!(basis, 1)
    basis.nprocessed += 1
    basis.nnonredundant += 1
    basis.nonredundant[basis.nnonredundant] = basis.nnonredundant
    basis.monoms[basis.nprocessed] = rexps
    basis.coeffs[basis.nprocessed] = rcoeffs
end

function divides_staircase(monom, staircase, ht)
    for m in staircase
        # if debug()
        #     # @warn "uwu"
        #     println("is $(ht.monoms[m]) divisible by $(ht.monoms[monom])")
        # end
        if is_monom_divisible(monom, m, ht)
            return true
        end
    end
    false
end

#------------------------------------------------------------------------------

function fglm_f4!(
    ring::PolyRing,
    basis::Basis{C},
    ht::MonomialHashtable,
    ord::AbstractMonomialOrdering
) where {C <: Coeff}
    newbasis = initialize_basis(ring, basis.nfilled, C)
    nextmonoms = initialize_nextmonomials(ht, ord)
    matrix = initialize_double_matrix(basis)
    staircase = MonomIdx[]

    while !isempty(nextmonoms)
        monom = nextmonomial!(nextmonoms)

        # if debug()
        #     # @warn "new iteration" ht.monoms[monom]
        #     # @warn "" ht
        #     println("#############")
        #     println(nextmonoms)
        #     println(newbasis.monoms, " ", newbasis.coeffs)
        # end

        if divides_staircase(monom, staircase, ht)
            continue
        end

        tobereduced = initialize_basis(ring, [[monom]], [C[1]])
        tobereduced.nfilled = 1

        # compute normal form
        f4_normalform!(ring, basis, tobereduced, ht)

        # if debug()
        #     println("Normal form ", tobereduced.monoms[1], " ", tobereduced.coeffs[1])
        # end

        # matrix left rows can express tobereduced?
        # reduces monom and tobereduced
        exists, relation = linear_relation!(ring, matrix, monom, tobereduced, ht)

        # if debug()
        #     # @warn "dump"
        #     dump(matrix, maxdepth=2)
        #     println(exists)
        # end

        # if linear relation between basis elements exists
        if exists
            lead = ht.monoms[monom]

            # if debug()
            #     # @error "produced element" lead
            # end

            add_generator!(newbasis, matrix, relation, ht, ord)
            push!(staircase, monom)
        else
            insertnexts!(nextmonoms, ht, monom)
        end
    end

    standardize_basis!(ring, newbasis, ht, ord)

    linbasis = extract_linear_basis(ring, matrix)

    newbasis, linbasis, ht
end

function fglm_f4(
    ring::PolyRing,
    basisexps::Vector{Vector{M}},
    basiscoeffs::Vector{Vector{C}},
    metainfo
) where {M, C <: Coeff}
    basis, pairset, ht = initialize_structs(ring, basisexps, basiscoeffs, metainfo)

    basis, linbasis, ht = fglm_f4!(ring, basis, ht, metainfo.computation_ord)

    export_basis_data(basis, ht)
end

#------------------------------------------------------------------------------

function extract_linear_basis(ring, matrix::DoubleMacaulayMatrix{C}) where {C}
    exps = Vector{Vector{MonomIdx}}(undef, matrix.nrrows)
    coeffs = Vector{Vector{C}}(undef, matrix.nrrows)

    for i in 1:(matrix.nrrows)
        exps[i] = matrix.rightrows[i]
        coeffs[i] = matrix.rightcoeffs[i]
        for j in 1:length(exps[i])
            exps[i][j] = matrix.rightcol2hash[exps[i][j]]
        end
    end

    linbasis = initialize_basis(ring, exps, coeffs)

    linbasis.nprocessed = length(exps)
    linbasis.nnonredundant = length(exps)
    linbasis.nonredundant = collect(1:(linbasis.nprocessed))

    linbasis
end

function _kbase(polynomials, kws)
    representation = select_polynomial_representation(polynomials, kws)
    ring, monoms, coeffs = convert_to_internal(representation, polynomials, kws)
    params = AlgorithmParameters(ring, kws)
    if isempty(monoms)
        @log level = -2 "Input consisting of zero polynomials."
        throw(DomainError("Input consisting of zero polynomials to Groebner.kbase."))
        return convert_to_output(ring, polynomials, monoms, coeffs, params)
    end
    if kws.check
        @log level = -2 "Checking if a Grobner basis"
        if !isgroebner(polynomials)
            throw(DomainError("Input is not a Groebner basis."))
        end
    end
    m, c = kbase_f4(ring, monoms, coeffs, params, representation)
    convert_to_output(ring, polynomials, m, c, params)
end

function kbase_f4(
    ring::PolyRing,
    basisexps::Vector{Vector{M}},
    basiscoeffs::Vector{Vector{C}},
    metainfo,
    representation
) where {M, C <: Coeff}
    # tablesize = select_tablesize(ring, basisexps)

    basis, pairset, ht = initialize_structs(ring, basisexps, basiscoeffs, metainfo)

    basis, linbasis, ht = fglm_f4!(ring, basis, ht, metainfo.computation_ord)

    export_basis_data(linbasis, ht)

    # #= set the logger =#
    # # prev_logger = Logging.global_logger(ConsoleLogger(stderr, loglevel))

    # # check && _check_isgroebner(basis)

    # #= extract ring information, exponents and coefficients
    #    from input polynomials =#
    # # Copies input, so that polynomials would not be changed itself.
    # ring, exps, coeffs =
    #     convert_to_internal(representation, basis, InputOrdering())

    # # metainfo = set_metaparameters(ring, InputOrdering(), false, :exact, rng)

    # iszerobasis = remove_zeros_from_input!(ring, exps, coeffs)
    # iszerobasis &&
    #     (throw(DomainError(basis, "Groebner.kbase does not work with such ideals, sorry")))

    # bexps, bcoeffs = kbase_f4(ring, exps, coeffs, metainfo)

    # # ordering in bexps here matches target ordering in metainfo

    # #= revert logger =#
    # # Logging.global_logger(prev_logger)

    # # ring contains ordering of computation, it is the requested ordering
    # #= convert result back to representation of input =#
    # convert_to_output(ring, basis, bexps, bcoeffs, metainfo)
end
