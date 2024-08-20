# This file is a part of Groebner.jl. License is GNU GPL v2.

# Groebner.jl interface goes in four levels of access:
# high level (0), intermediate level (1), low level (2), F4 level (3).
#
# Level (0) works with high-level polynomials from some frontend.
# Level (1) works with exponent vectors and coefficients.
# Level (2) works with internal representations of monomials and coefficients.
# Level (3) works with internal F4 data structures (matrices, hashtables, etc).
#
# User interface is levels (0) and (1). Useful work is done on levels (2) and
# (3). F4 runs on level (3). A series of conversions carries data across levels.
# Generally, input and output of a function must live on the same level.

###
# Backend for `groebner`

# polynomials => polynomials
function groebner0(polynomials::AbstractVector, options::KeywordArguments)
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    gb_monoms, gb_coeffs = _groebner1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    result
end

# (exponent vectors, coefficients) => (exponent vectors, coefficients)
function groebner1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    ring, monoms, coeffs = ir_ensure_assumptions(ring, monoms, coeffs)
    _groebner1(ring, monoms, coeffs, options)
end

function _groebner1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    options::KeywordArguments
) where {I <: Integer, C <: Coeff}
    try
        params = AlgorithmParameters(ring, options)
        return __groebner1(ring, monoms, coeffs, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent.""" maxlog = 1
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __groebner1(ring, monoms, coeffs, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __groebner1(
    ring::PolyRing,
    monoms::Vector{Vector{Vector{I}}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {I <: Integer, C <: Coeff}
    @invariant ir_is_valid(ring, monoms, coeffs)
    _, ring2, monoms2, coeffs2 =
        io_convert_ir_to_internal(ring, monoms, coeffs, params, params.representation)
    gb_monoms2, gb_coeffs2 = groebner2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = io_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)
    gb_monoms, gb_coeffs
end

# internal structs => internal structs
function groebner2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: Coeff}
    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return [monoms[1]], [coeffs[1]]
    end
    monoms, coeffs = _monoms, _coeffs

    if params.homogenize
        _, ring, monoms, coeffs = homogenize_generators!(ring, monoms, coeffs, params)
    end

    gb_monoms, gb_coeffs = _groebner2(ring, monoms, coeffs, params)

    if params.homogenize
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end

    gb_monoms, gb_coeffs
end

###
# Groebner basis over Z_p. Calls F4 directly.

function _groebner2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffZp}
    @log :misc "Backend: F4 over Z_$(ring.ch)"
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    f4!(ring, basis, pairset, hashtable, params)
    gbmonoms, gbcoeffs = basis_export_data(basis, hashtable)
    gbmonoms, gbcoeffs
end

###
# Groebner basis over Q.

# GB over the rationals uses modular computation.
function _groebner2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    if params.modular_strategy === :learn_and_apply
        if (
            params.threaded_multimodular === :yes || params.threaded_multimodular === :auto
        ) && nthreads() > 1
            _groebner_learn_and_apply_threaded(ring, monoms, coeffs, params)
        else
            _groebner_learn_and_apply(ring, monoms, coeffs, params)
        end
    else
        @assert params.modular_strategy === :classic_modular
        _groebner_classic_modular(ring, monoms, coeffs, params)
    end
end

###
# Learn & Apply startegy

# The next batchsize is a multiple of the previous one aligned to some nice
# power of two.
function get_next_batchsize(
    primes_used::Int,
    prev_batchsize::Int,
    batchsize_scaling::Float64
)
    new_batchsize = if prev_batchsize == 1
        4
    else
        max(4, round(Int, primes_used * batchsize_scaling))
    end
    if new_batchsize > 2
        new_batchsize::Int = align_up(new_batchsize, 4)
    end
    max(new_batchsize, prev_batchsize)
end

function _groebner_learn_and_apply(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}

    # Initialize supporting structs
    state = GroebnerState{BigInt, C, Int32}(params)
    basis, pairset, hashtable, permutation =
        f4_initialize_structs(ring, monoms, coeffs, params, make_monic=false)

    basis_zz = clear_denominators!(state.buffer, basis, deepcopy=false)

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_lucky_prime!(luckyprimes)

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)

    params_zp = params_mod_p(params, prime)
    trace = trace_initialize(
        ring_ff,
        basis_deepcopy(basis_ff),
        basis_ff,
        hashtable,
        permutation,
        params_zp
    )

    f4_learn!(trace, ring_ff, trace.gb_basis, pairset, hashtable, params_zp)

    # TODO: no need to deepcopy!
    push!(state.gb_coeffs_ff_all, deepcopy(trace.gb_basis.coeffs))

    # Reconstruct coefficients and write results to the accumulator.
    full_simultaneous_crt_reconstruct!(state, luckyprimes)

    success_reconstruct = full_rational_reconstruct!(state, luckyprimes, params.use_flint)

    correct_basis = false
    if success_reconstruct
        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            trace.gb_basis,
            hashtable,
            params
        )
        # At this point, if the constructed basis is correct, we return it.
        if correct_basis
            gb_monoms, _ = basis_export_data(trace.gb_basis, hashtable)
            gb_coeffs_qq = state.gb_coeffs_qq
            return gb_monoms, gb_coeffs_qq
        end
    end

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    primes_used = 1
    batchsize = 1
    batchsize_scaling = 0.10

    # CRT and rational reconstrction settings
    tot, indices_selection = gb_modular_select_indices0(state.gb_coeffs_zz, params)

    # Initialize partial CRT reconstruction
    partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

    # Initialize a tracer that computes the bases in batches of 4
    trace_4x = trace_copy(trace, CompositeNumber{4, Int32})

    iters = 0
    while !correct_basis
        if iszero(batchsize % 4) && params.batched
            for j in 1:4:batchsize
                prime_4x = ntuple(_ -> Int32(next_lucky_prime!(luckyprimes)), 4)

                # Perform reduction modulo primes and store result in basis_ff_4x
                ring_ff_4x, basis_ff_4x =
                    reduce_modulo_p_in_batch!(state.buffer, ring, basis_zz, prime_4x)
                params_zp_4x = params_mod_p(
                    params,
                    CompositeNumber{4, Int32}(prime_4x),
                    using_wide_type_for_coeffs=false
                )
                trace_4x.buf_basis = basis_ff_4x
                trace_4x.ring = ring_ff_4x

                f4_apply!(trace_4x, ring_ff_4x, trace_4x.buf_basis, params_zp_4x)
                gb_coeffs_1, gb_coeffs_2, gb_coeffs_3, gb_coeffs_4 =
                    io_unpack_composite_coefficients(trace_4x.gb_basis.coeffs)

                # TODO: (I) This causes unnecessary conversions of arrays.
                push!(state.gb_coeffs_ff_all, gb_coeffs_1)
                push!(state.gb_coeffs_ff_all, gb_coeffs_2)
                push!(state.gb_coeffs_ff_all, gb_coeffs_3)
                push!(state.gb_coeffs_ff_all, gb_coeffs_4)
                primes_used += 4
            end
        else
            for j in 1:batchsize
                prime = next_lucky_prime!(luckyprimes)

                # Perform reduction modulo prime and store result in basis_ff
                ring_ff, basis_ff =
                    reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
                params_zp = params_mod_p(params, prime)

                trace.buf_basis = basis_ff
                trace.ring = ring_ff

                f4_apply!(trace, ring_ff, trace.buf_basis, params_zp)

                push!(state.gb_coeffs_ff_all, deepcopy(trace.gb_basis.coeffs))

                if !majority_vote!(state, trace.gb_basis, params)
                    continue
                end
                primes_used += 1
            end
        end

        partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

        success_reconstruct = partial_rational_reconstruct!(
            state,
            luckyprimes,
            indices_selection,
            params.use_flint
        )

        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                iters += 1
                batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
                continue
            end
        end

        # Perform full reconstruction
        full_simultaneous_crt_reconstruct!(state, luckyprimes)

        success_reconstruct =
            full_rational_reconstruct!(state, luckyprimes, params.use_flint)

        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            trace.gb_basis,
            hashtable,
            params
        )

        iters += 1
        batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
    end

    # Construct the output basis.
    # Take monomials from the basis modulo a prime
    gb_monoms, _ = basis_export_data(trace.gb_basis, hashtable)
    # Take coefficients from the reconstructed basis
    gb_coeffs_qq = state.gb_coeffs_qq

    return gb_monoms, gb_coeffs_qq
end

###
# Threaded Learn & Apply startegy

function _groebner_learn_and_apply_threaded(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    if nthreads() == 1
        @log :info "Using threaded backend with nthreads() == 1, how did we end up here?"
    end

    # Initialize supporting structs
    state = GroebnerState{BigInt, C, Int32}(params)
    basis, pairset, hashtable, permutation =
        f4_initialize_structs(ring, monoms, coeffs, params, make_monic=false)

    basis_zz = clear_denominators!(state.buffer, basis, deepcopy=false)

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_lucky_prime!(luckyprimes)

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)

    params_zp = params_mod_p(params, prime)
    trace = trace_initialize(
        ring_ff,
        basis_deepcopy(basis_ff),
        basis_ff,
        hashtable,
        permutation,
        params_zp
    )

    f4_learn!(trace, ring_ff, trace.gb_basis, pairset, hashtable, params_zp)

    # TODO: no need to deepcopy!
    push!(state.gb_coeffs_ff_all, deepcopy(trace.gb_basis.coeffs))

    # Reconstruct coefficients and write results to the accumulator.
    full_simultaneous_crt_reconstruct!(state, luckyprimes)

    success_reconstruct = full_rational_reconstruct!(state, luckyprimes, params.use_flint)

    correct_basis = false
    if success_reconstruct
        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            trace.gb_basis,
            hashtable,
            params
        )
        if correct_basis
            gb_monoms, _ = basis_export_data(trace.gb_basis, hashtable)
            gb_coeffs_qq = state.gb_coeffs_qq
            return gb_monoms, gb_coeffs_qq
        end
    end

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    primes_used = 1
    batchsize = align_up(min(32, 4 * nthreads()), 4)
    batchsize_scaling = 0.10

    # CRT and rational reconstrction settings
    tot, indices_selection = gb_modular_select_indices0(state.gb_coeffs_zz, params)

    # Initialize partial CRT reconstruction
    partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

    # Initialize a tracer that computes the bases in batches of 4
    trace_4x = trace_copy(trace, CompositeNumber{4, Int32})

    # Thread buffers
    threadbuf_trace_4x = map(_ -> trace_deepcopy(trace_4x), 1:nthreads())
    threadbuf_gb_coeffs =
        Vector{Vector{Tuple{Int32, Vector{Vector{Int32}}}}}(undef, nthreads())
    for i in 1:nthreads()
        threadbuf_gb_coeffs[i] = Vector{Tuple{Int, Vector{Vector{Int32}}}}()
    end
    threadbuf_bigint_buffer = map(_ -> CoefficientBuffer(), 1:nthreads())
    threadbuf_params = map(_ -> deepcopy(params), 1:nthreads())

    iters = 0
    while !correct_basis
        @invariant iszero(batchsize % 4)

        threadbuf_primes = map(_ -> Int32(next_lucky_prime!(luckyprimes)), 1:batchsize)
        for i in 1:nthreads()
            empty!(threadbuf_gb_coeffs[i])
        end

        Base.Threads.@threads :static for j in 1:4:batchsize
            t_id = threadid()
            threadlocal_trace_4x = threadbuf_trace_4x[t_id]
            threadlocal_prime_4x = ntuple(k -> threadbuf_primes[j + k - 1], 4)
            threadlocal_bigint_buffer = threadbuf_bigint_buffer[t_id]
            threadlocal_params = threadbuf_params[t_id]

            ring_ff_4x, basis_ff_4x = reduce_modulo_p_in_batch!(
                threadlocal_bigint_buffer,  # is modified
                ring,                       # is not modified
                basis_zz,                   # is not modified
                threadlocal_prime_4x        # is not modified
            )
            threadlocal_params_zp_4x = params_mod_p(
                threadlocal_params,               # can be modified later
                CompositeNumber{4, Int32}(threadlocal_prime_4x),
                using_wide_type_for_coeffs=false
            )
            threadlocal_trace_4x.buf_basis = basis_ff_4x
            threadlocal_trace_4x.ring = ring_ff_4x

            f4_apply!(
                threadlocal_trace_4x,
                ring_ff_4x,
                threadlocal_trace_4x.buf_basis,
                threadlocal_params_zp_4x
            )

            gb_coeffs_1, gb_coeffs_2, gb_coeffs_3, gb_coeffs_4 =
                io_unpack_composite_coefficients(threadlocal_trace_4x.gb_basis.coeffs)

            push!(threadbuf_gb_coeffs[t_id], (threadlocal_prime_4x[1], gb_coeffs_1))
            push!(threadbuf_gb_coeffs[t_id], (threadlocal_prime_4x[2], gb_coeffs_2))
            push!(threadbuf_gb_coeffs[t_id], (threadlocal_prime_4x[3], gb_coeffs_3))
            push!(threadbuf_gb_coeffs[t_id], (threadlocal_prime_4x[4], gb_coeffs_4))
        end
        primes_used += batchsize

        threadbuf_gb_coeffs_union = reduce(vcat, threadbuf_gb_coeffs)

        sort!(threadbuf_gb_coeffs_union, by=first, rev=true)
        for (_, coeffs_ff_) in threadbuf_gb_coeffs_union
            push!(state.gb_coeffs_ff_all, coeffs_ff_)
        end

        partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

        success_reconstruct = partial_rational_reconstruct!(
            state,
            luckyprimes,
            indices_selection,
            params.use_flint
        )

        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                iters += 1
                batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
                continue
            end
        end

        # Perform full reconstruction
        full_simultaneous_crt_reconstruct!(state, luckyprimes)
        success_reconstruct =
            full_rational_reconstruct!(state, luckyprimes, params.use_flint)

        # This should happen rarely
        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            trace.gb_basis,
            hashtable,
            params
        )

        iters += 1
        batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
    end

    gb_monoms, _ = basis_export_data(trace.gb_basis, hashtable)
    gb_coeffs_qq = state.gb_coeffs_qq

    return gb_monoms, gb_coeffs_qq
end

###
# Classic multi-modular strategy

function _groebner_classic_modular(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}

    # Initialize supporting structs
    state = GroebnerState{BigInt, C, CoeffModular}(params)
    basis, pairset, hashtable =
        f4_initialize_structs(ring, monoms, coeffs, params, make_monic=false)

    basis_zz = clear_denominators!(state.buffer, basis, deepcopy=false)

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_lucky_prime!(luckyprimes)

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)

    params_zp = params_mod_p(params, prime)
    f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)
    # NOTE: basis_ff may not own its coefficients, one should not mutate it
    # directly further in the code

    push!(state.gb_coeffs_ff_all, basis_ff.coeffs)

    # Reconstruct coefficients and write results to the accumulator.
    full_simultaneous_crt_reconstruct!(state, luckyprimes)

    success_reconstruct = full_rational_reconstruct!(state, luckyprimes, params.use_flint)

    correct_basis = false
    if success_reconstruct
        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            basis_ff,
            hashtable,
            params
        )
        if correct_basis
            gb_monoms, _ = basis_export_data(basis_ff, hashtable)
            gb_coeffs_qq = state.gb_coeffs_qq
            return gb_monoms, gb_coeffs_qq
        end
    end

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    primes_used = 1
    batchsize = 1
    batchsize_scaling = 0.10

    # CRT and rational reconstrction settings
    tot, indices_selection = gb_modular_select_indices0(state.gb_coeffs_zz, params)

    partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

    iters = 0
    while !correct_basis
        for j in 1:batchsize
            prime = next_lucky_prime!(luckyprimes)

            # Perform reduction modulo prime and store result in basis_ff
            ring_ff, basis_ff =
                reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
            params_zp = params_mod_p(params, prime)

            f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)

            push!(state.gb_coeffs_ff_all, basis_ff.coeffs)

            if !majority_vote!(state, basis_ff, params)
                continue
            end
            primes_used += 1
        end

        partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

        success_reconstruct = partial_rational_reconstruct!(
            state,
            luckyprimes,
            indices_selection,
            params.use_flint
        )

        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                iters += 1
                batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
                continue
            end
        end

        # Perform full reconstruction
        full_simultaneous_crt_reconstruct!(state, luckyprimes)
        success_reconstruct =
            full_rational_reconstruct!(state, luckyprimes, params.use_flint)

        if !success_reconstruct
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        correct_basis = correctness_check!(
            state,
            luckyprimes,
            ring_ff,
            basis,
            basis_zz,
            basis_ff,
            hashtable,
            params
        )

        iters += 1
        batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
    end

    gb_monoms, _ = basis_export_data(basis_ff, hashtable)
    gb_coeffs_qq = state.gb_coeffs_qq

    return gb_monoms, gb_coeffs_qq
end
