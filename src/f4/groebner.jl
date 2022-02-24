
#######################################

# Stores parameters for groebner algorithm
struct Metainfo
    # if set, then use fglm algorithm for order conversion
    usefglm::Bool
    # output polynomials order
    targetord::Symbol
    # order for computation
    computeord::Symbol
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

    @info "Computing in $computeord order, result is in $targetord order"
    @info "Using fglm: $usefglm"

    Metainfo(usefglm, targetord, computeord)
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
            rng::Rng) where {Rng<:Random.AbstractRNG}

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
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{Rational{BigInt}}},
            reduced::Bool,
            randomized::Bool,
            rng::Rng) where {Rng<:Random.AbstractRNG}

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


    gbcoeffs_accum = Vector{Vector{BigInt}}(undef, 0)
    gbcoeffs_qq = Vector{Vector{Rational{BigInt}}}(undef, 0)

    init_coeffs_zz = scale_denominators(coeffs)

    # @warn "scaled" coeffs init_coeffs_zz

    modulo = BigInt(1)
    prime = nextluckyprime(init_coeffs_zz)
    moduli = Int[prime]
    goodprime = nextgoodprime(init_coeffs_zz, moduli)

    # @error "initial primes" prime goodprime

    ring_ff, init_coeffs_ff = reduce_modulo(init_coeffs_zz, ring, prime)

    # @warn "reduced" ring_ff init_coeffs_ff init_coeffs_zz

    # TODO: 2^16

    init_gens_temp_ff, ht = initialize_structures(ring_ff, exps,
                                                    init_coeffs_zz, init_coeffs_ff,
                                                    rng, tablesize)
    gens_ff = copy_basis(init_gens_temp_ff)

    @assert gens_ff.ch == prime == ring_ff.ch
    # @warn "after init" gens_ff

    i = 1
    while true
        @info "$i: selected lucky prime $prime"

        @assert ring_ff.ch == prime
        @assert gens_ff.ch == prime
        @assert prime < 2^32
        # @error "Asserts ok!"

        # @info "initial gens" gens_ff.coeffs
        #@error "passing gens to f4"
        #println(gens_ff.lead)
        #println("OK!")

        f4!(ring_ff, gens_ff, ht, reduced)

        #@info "basis computed" gens_ff.coeffs

        @info "CRT modulo ($modulo, $(ring_ff.ch))"
        reconstruct_crt!(gbcoeffs_accum, modulo, gens_ff.coeffs, ring_ff)

        #@info "reconstructed #1" gbcoeffs_accum

        @info "Reconstructing modulo $modulo"
        reconstruct_modulo!(gbcoeffs_qq, gbcoeffs_accum, modulo)

        #@info "reconstructed #2" gbcoeffs_qq
        #println(ht.exponents[1:10])

        #@error "passing gens to correctness check"
        #println(map(bitstring, gens_ff.lead))
        #println([ht.exponents[i[1]] for i in gens_ff.gens])
        #println("OK!")

        buf_ff = copy_basis(init_gens_temp_ff)
        if correctness_check!(init_coeffs_zz, ring_ff, gens_ff,
                                        buf_ff, ht, gbcoeffs_qq, gbcoeffs_accum,
                                        modulo, randomized, goodprime)
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
    gbexps = hash_to_exponents(gens_ff, ht)
    gbexps, gbcoeffs_qq
end

#######################################
# Finite field isgroebner

# UWU!
function isgroebner_ff(
            ring::PolyRing,
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{UInt64}},
            rng)

    isgroebner_f4(ring, exps, coeffs, rng)
end

#######################################
# Rational field groebner

# TODO
function isgroebner_qq(
            ring::PolyRing,
            exps::Vector{Vector{Vector{UInt16}}},
            coeffs::Vector{Vector{Rational{BigInt}}},
            randomized::Bool,
            rng)

    if randomized
        coeffs_zz = scale_denominators(coeffs)
        goodprime = nextgoodprime(coeffs_zz, Int[], 2^30 + 3)
        ring_ff, coeffs_ff = reduce_modulo(coeffs_zz, ring, goodprime)
        isgroebner_f4(ring_ff, exps, coeffs_ff, rng)
    else
        error("Sorry, not randomized version is not implemented yet.")
    end
end
