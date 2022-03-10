
#######################################

# Stores parameters for groebner algorithm
struct GroebnerMetainfo
    # if set, then use fglm algorithm for order conversion
    usefglm::Bool
    # output polynomials order
    targetord::Symbol
    # order for computation
    computeord::Symbol

    # correctness checks levels
    heuristiccheck::Bool
    randomizedcheck::Bool
    guaranteedcheck::Bool
end

function set_metaparameters(ring, ordering, certify, forsolve)
    usefglm = false
    targetord = :lex
    computeord = :lex

    if forsolve
        targetord = :lex
        usefglm = true
        if ordering in (:deglex, :degrevlex)
            computeord = ordering
        else
            computeord = :degrevlex
        end
        # TODO:
        computeord = :lex
    else
        if ordering == :input
            ordering = ring.ord
        end
        targetord = ordering
        usefglm = false
        computeord = targetord
    end

    heuristiccheck = true
    randomizedcheck = true
    if certify
        guaranteedcheck = true
    else
        guaranteedcheck = false
    end

    @info "Computing in $computeord order, result is in $targetord order"
    @info "Using fglm: $usefglm"

    GroebnerMetainfo(usefglm, targetord, computeord,
                        heuristiccheck, randomizedcheck, guaranteedcheck)
end

function assure_ordering!(ring, exps, coeffs, metainfo)
    if ring.ord != metainfo.computeord
        sort_input_to_change_ordering(exps, coeffs, metainfo.computeord)
    end
    ring.ord = metainfo.computeord
end

# Mutate everything!

#######################################
# Finite field groebner

function groebner_ff(
            ring::PolyRing,
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{UInt64}},
            reduced::Bool,
            rng::Rng,
            meta::GroebnerMetainfo) where {Rng<:Random.AbstractRNG}

    # TODO
    tablesize = 2^8
    if ring.nvars > 3
        tablesize = 2^10
    end
    if ring.nvars > 5
        tablesize = 2^14
    end
    tablesize = 2^16

    # specialize on ordering (not yet)
    # groebner_ff(ring, exps, coeffs, reduced, rng, Val(ring.ord))
    basis, ht = initialize_structures(
                        ring, exps, coeffs, rng, tablesize)

    f4!(ring, basis, ht, reduced)

    gbexps = hash_to_exponents(basis, ht)
    gbexps, basis.coeffs
end

#######################################
# Rational field groebner

function groebner_qq(
            ring::PolyRing,
            exps::Vector{Vector{ExponentVector}},
            coeffs::Vector{Vector{CoeffQQ}},
            reduced::Bool,
            rng::Rng,
            meta::GroebnerMetainfo) where {Rng<:Random.AbstractRNG}

    # we can mutate coeffs and exps here
    #=
        TODO
    =#
    tablesize = 2^8
    if ring.nvars > 3
        tablesize = 2^10
    end
    if ring.nvars > 5
        tablesize = 2^14
    end
    tablesize = 2^16

    # to store integer and rational coefficients of groebner basis
    gb_coeffs_accum = Vector{Vector{CoeffZZ}}(undef, 0)
    gb_coeffs_qq = Vector{Vector{CoeffQQ}}(undef, 0)

    # to store integer coefficients of input
    init_coeffs_zz = scale_denominators(coeffs)

    # selecting first lucky prime and first good prime
    modulo = BigInt(1)
    prime = nextluckyprime(init_coeffs_zz)
    moduli = Int[prime]
    goodprime = nextgoodprime(init_coeffs_zz, moduli)

    # first reduction
    ring_ff, init_coeffs_ff = reduce_modulo(init_coeffs_zz, ring, prime)

    # TODO: 2^16

    # temporary basis for initial generators in finite field
    init_gens_temp_ff, ht = initialize_structures(ring_ff, exps,
                                            coeffs, init_coeffs_zz,
                                            init_coeffs_ff, rng, tablesize)
    gens_ff = copy_basis(init_gens_temp_ff)

    @assert gens_ff.ch == prime == ring_ff.ch

    i = 1
    while true
        @info "$i: selected lucky prime $prime"

        @assert ring_ff.ch == prime
        @assert gens_ff.ch == prime
        @assert prime < 2^32

        # compute groebner basis in finite field
        f4!(ring_ff, gens_ff, ht, reduced)

        # reconstruct into integers
        @info "CRT modulo ($modulo, $(ring_ff.ch))"
        reconstruct_crt!(gb_coeffs_accum, modulo, gens_ff.coeffs, ring_ff)

        # reconstruct into rationals
        @info "Reconstructing modulo $modulo"
        reconstruct_modulo!(gb_coeffs_qq, gb_coeffs_accum, modulo)

        buf_ff = copy_basis(init_gens_temp_ff)
        if correctness_check!(init_gens_temp_ff, coeffs, init_coeffs_zz,
                                ring_ff, gens_ff, buf_ff, ht, gb_coeffs_qq,
                                gb_coeffs_accum, modulo, goodprime, meta,
                                tablesize, rng)
            @info "Success!"
            break
        end

        prime = nextluckyprime(init_coeffs_zz, prime)
        push!(moduli, prime)
        goodprime = nextgoodprime(init_coeffs_zz, moduli, goodprime)

        reduce_modulo!(init_coeffs_zz, prime, ring_ff, init_gens_temp_ff)
        gens_ff = copy_basis(init_gens_temp_ff)
        normalize_basis!(gens_ff)

        i += 1
    end

    # normalize_coeffs!(gbcoeffs_qq)
    gb_exps = hash_to_exponents(gens_ff, ht)
    gb_exps, gb_coeffs_qq
end

#######################################
# Finite field isgroebner

# UWU!
function isgroebner_ff(
            ring::PolyRing,
            exps::Vector{Vector{ExponentVector}},
            coeffs::Vector{Vector{CoeffFF}},
            rng,
            meta)

    isgroebner_f4(ring, exps, coeffs, rng)
end

#######################################
# Rational field groebner

# TODO
function isgroebner_qq(
            ring::PolyRing,
            exps::Vector{Vector{ExponentVector}},
            coeffs::Vector{Vector{CoeffQQ}},
            rng,
            meta)

    # if randomized result is ok
    if !meta.guaranteedcheck
        coeffs_zz = scale_denominators(coeffs)
        goodprime = nextgoodprime(coeffs_zz, Int[], 2^30 + 3)
        ring_ff, coeffs_ff = reduce_modulo(coeffs_zz, ring, goodprime)
        isgroebner_f4(ring_ff, exps, coeffs_ff, rng)
    else # if proved result is needed, compute in rationals
        isgroebner_f4(ring, exps, coeffs, rng)
    end
end
