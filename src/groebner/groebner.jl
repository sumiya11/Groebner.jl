# Backend for `groebner`

# Proxy function for handling exceptions.
function _groebner(polynomials, kws::KeywordsHandler)
    # We try to select an efficient internal polynomial representation, i.e., a
    # suitable representation of monomials and coefficients.
    polynomial_repr = select_polynomial_representation(polynomials, kws)
    try
        # The backend is wrapped in a try/catch to catch exceptions that one can
        # hope to recover from (and, perhaps, restart the computation with safer
        # parameters).
        return _groebner(polynomials, kws, polynomial_repr)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log level = 1 """
            Possible overflow of exponent vector detected. 
            Restarting with at least $(32) bits per exponent."""
            polynomial_repr =
                select_polynomial_representation(polynomials, kws, hint=:large_exponents)
            return _groebner(polynomials, kws, polynomial_repr)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function _groebner(polynomials, kws::KeywordsHandler, representation)
    # Extract ring information, exponents, and coefficients from input
    # polynomials. Convert these to an internal polynomial representation. 
    # NOTE: This must copy the input, so that input `polynomials` is never
    # modified.
    # NOTE: The body of this function is type-unstable (by design)
    ring, var_to_index, monoms, coeffs =
        convert_to_internal(representation, polynomials, kws)

    # Check and set parameters and monomial ordering
    params = AlgorithmParameters(ring, representation, kws)
    ring, _ = set_monomial_ordering!(ring, var_to_index, monoms, coeffs, params)

    # Fast path for the input of zeros
    if isempty(monoms)
        @log level = -2 "Input consisting of zero polynomials. Returning zero."
        return convert_to_output(ring, polynomials, monoms, coeffs, params)
    end

    if params.homogenize
        # this also performs saturation w.r.t. the homogenizing variable
        _, ring, monoms, coeffs = homogenize_generators!(ring, monoms, coeffs, params)
    end

    # Compute a groebner basis!
    gbmonoms, gbcoeffs = _groebner(ring, monoms, coeffs, params)

    if params.homogenize
        ring, gbmonoms, gbcoeffs =
            dehomogenize_generators!(ring, gbmonoms, gbcoeffs, params)
    end

    # Convert result back to the representation of input
    basis = convert_to_output(ring, polynomials, gbmonoms, gbcoeffs, params)

    print_performance_counters(params.statistics)
    print_statistics(params.statistics)

    basis
end

###
# Groebner basis over Z_p.
# Just calls f4 directly.

@timeit function _groebner(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffFF}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    @log level = -1 "Backend: F4 over Z_$(ring.ch)"
    # NOTE: the sorting of input polynomials is not deterministic across
    # different Julia versions when sorting only w.r.t. the leading term
    basis, pairset, hashtable = initialize_structs(ring, monoms, coeffs, params)
    tracer = TinyTraceF4()
    f4!(ring, basis, pairset, hashtable, tracer, params)
    # Extract monomials and coefficients from basis and hashtable
    gbmonoms, gbcoeffs = export_basis_data(basis, hashtable)
    gbmonoms, gbcoeffs
end

###
# Groebner basis over Q.

# GB over the rationals uses modular computation.
@timeit function _groebner(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    if params.modular_strategy === :learn_and_apply
        _groebner_learn_and_apply(ring, monoms, coeffs, params)
    else
        @assert params.modular_strategy === :classic_modular
        _groebner_classic_modular(ring, monoms, coeffs, params)
    end
end

function get_next_batchsize(batchsize, batchsize_multiplier)
    max(batchsize + 1, round(Int, batchsize * batchsize_multiplier))
end

function _groebner_learn_and_apply(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    @log level = -1 "Backend: learn & apply multi-modular F4"

    # Initialize supporting structs
    state = GroebnerState{BigInt, C, CoeffModular}(params)
    # Initialize F4 structs
    basis, pairset, hashtable, permutation =
        initialize_structs(ring, monoms, coeffs, params, normalize_input=false)
    crt_algorithm = params.crt_algorithm

    # Scale the input coefficients to integers to speed up the subsequent search
    # for lucky primes
    @log level = -5 "Input polynomials" basis
    @log level = -2 "Clearing the denominators of the input polynomials"
    basis_zz = clear_denominators!(state.buffer, basis, deepcopy=false)
    @log level = -5 "Integer coefficients are" basis_zz.coeffs

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_lucky_prime!(luckyprimes)
    @log level = -2 "The first lucky prime is $prime"
    @log level = -2 "Reducing input generators modulo $prime"

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
    @log level = -5 "Reduced coefficients are" basis_ff.coeffs

    @log level = -5 "Before F4" basis_ff
    params_zp = params_mod_p(params, prime)
    trace = initialize_trace_f4(
        ring_ff,
        deepcopy_basis(basis_ff),
        basis_ff,
        hashtable,
        permutation,
        params_zp
    )

    f4_learn!(trace, ring_ff, trace.gb_basis, pairset, hashtable, params_zp)
    @log level = -5 "After F4:" trace.gb_basis

    push!(state.gb_coeffs_ff_all, deepcopy(trace.gb_basis.coeffs))

    # Reconstruct coefficients and write results to the accumulator.
    # CRT reconstrction is trivial here.
    @log level = -2 "Reconstructing coefficients from Z_$prime to QQ"
    if crt_algorithm === :incremental
        full_incremental_crt_reconstruct!(state, luckyprimes)
    else
        full_simultaneous_crt_reconstruct!(state, luckyprimes)
    end
    success_reconstruct = full_rational_reconstruct!(state, luckyprimes)
    @log level = -5 "Reconstructed coefficients" state.gb_coeffs_qq
    @log level = -2 "Successfull reconstruction: $success_reconstruct"

    correct_basis = false
    if success_reconstruct
        @log level = -2 "Verifying the correctness of reconstruction"
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
        @log level = -2 "Passed correctness check: $correct_basis"
        # At this point, if the constructed basis is correct, we return it.
        if correct_basis
            # take monomials from the basis modulo a prime
            gb_monoms, _ = export_basis_data(trace.gb_basis, hashtable)
            # take coefficients from the reconstrcted basis
            gb_coeffs_qq = state.gb_coeffs_qq
            return gb_monoms, gb_coeffs_qq
        end
    end

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    batchsize = 1
    batchsize_multiplier = 1.4
    @log level = -2 """
      Preparing to compute bases in batches.. 
      The initial size of the batch is $batchsize. 
      The size increases in a geometric progression. 
      The batch size multiplier is $batchsize_multiplier.
      """

    # CRT and rational reconstrction settings
    indices_selection = Vector{Tuple{Int, Int}}(undef, length(state.gb_coeffs_zz))
    k = 1
    for i in 1:length(state.gb_coeffs_zz)
        l = length(state.gb_coeffs_zz[i])
        nl = max(isqrt(l) - 1, 1)
        if isone(l)
            continue
        end
        for j in 1:nl
            if k > length(indices_selection)
                resize!(indices_selection, 2 * length(indices_selection))
            end
            indices_selection[k] = (i, rand(2:l))
            k += 1
        end
    end
    resize!(indices_selection, k - 1)
    unique!(indices_selection)
    @log level = -2 "" length(indices_selection) length(state.gb_coeffs_zz) sum(
        map(length, state.gb_coeffs_zz)
    )

    if crt_algorithm === :incremental
        partial_incremental_crt_reconstruct!(state, luckyprimes, indices_selection)
    else
        partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)
    end

    iters = 0
    while !correct_basis
        @log level = -2 "Iteration # $iters of modular Groebner"

        for j in 1:batchsize
            prime = next_lucky_prime!(luckyprimes)
            @log level = -3 "The lucky prime is $prime"
            @log level = -3 "Reducing input generators modulo $prime"

            # Perform reduction modulo prime and store result in basis_ff
            ring_ff, basis_ff =
                reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
            params_zp = params_mod_p(params, prime)

            trace.buf_basis = basis_ff
            trace.ring = ring_ff

            f4_apply!(trace, ring_ff, trace.buf_basis, params_zp)

            push!(state.gb_coeffs_ff_all, deepcopy(trace.gb_basis.coeffs))

            if !majority_vote!(state, trace.gb_basis, nothing, params)
                @log level = -3 "Majority vote is not conclusive, aborting reconstruction!"
                continue
            end

            # @log level = -2 "Reconstructing coefficients using CRT"
            # @log level = -4 "Reconstructing coefficients from Z_$(luckyprimes.modulo) * Z_$(prime) to Z_$(luckyprimes.modulo * prime)"
            if crt_algorithm === :incremental
                partial_incremental_crt_reconstruct!(state, luckyprimes, indices_selection)
            end
        end
        # @log level = -2 "Reconstructing coefficients using CRT"
        if crt_algorithm === :simultaneous
            partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)
        end

        @log level = -2 "Reconstructing coefficients to QQ"
        @log level = -4 "Reconstructing coefficients from Z_$(luckyprimes.modulo * prime) to QQ"
        success_reconstruct =
            partial_rational_reconstruct!(state, luckyprimes, indices_selection)
        @log level = -2 "Reconstruction successfull: $success_reconstruct"
        @log level = -2 """
          Used $(length(luckyprimes.primes)) primes in total over $(iters + 1) iterations.
          The current batch size is $batchsize.
          """

        # @info "" state.gb_coeffs_ff_all
        # @info "" state.selected_coeffs_zz
        # @info "" state.selected_prev_coeffs_zz
        # @info "" state.selected_coeffs_qq

        if !success_reconstruct
            @log level = -2 "Partial rational reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                @log level = -2 "Heuristic check failed for partial reconstruction"
                iters += 1
                batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
                continue
            end
        end

        # Perform full reconstruction
        @log level = -2 "Performing full CRT and rational reconstruction.."
        if crt_algorithm === :incremental
            full_incremental_crt_reconstruct!(state, luckyprimes)
        else
            full_simultaneous_crt_reconstruct!(state, luckyprimes)
        end
        success_reconstruct = full_rational_reconstruct!(state, luckyprimes)

        if !success_reconstruct
            @log level = -2 "Full reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
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
        batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
    end

    @log level = -2 "Correctness check passed!"
    @log level = -2 "Used $(length(luckyprimes.primes)) primes in total over $(iters) iterations"

    # Construct the output basis.
    # Take monomials from the basis modulo a prime
    gb_monoms, _ = export_basis_data(trace.gb_basis, hashtable)
    # Take coefficients from the reconstructed basis
    gb_coeffs_qq = state.gb_coeffs_qq

    return gb_monoms, gb_coeffs_qq
end

function _groebner_classic_modular(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    @log level = -1 "Backend: classic multi-modular F4"

    # Initialize supporting structs
    state = GroebnerState{BigInt, C, CoeffModular}(params)
    # Initialize F4 structs
    basis, pairset, hashtable =
        initialize_structs(ring, monoms, coeffs, params, normalize_input=false)
    tracer = TinyTraceF4(params)
    crt_algorithm = params.crt_algorithm

    # Scale the input coefficients to integers to speed up the subsequent search
    # for lucky primes
    @log level = -5 "Input polynomials" basis
    @log level = -2 "Clearing the denominators of the input polynomials"
    basis_zz = clear_denominators!(state.buffer, basis, deepcopy=false)
    @log level = -5 "Integer coefficients are" basis_zz.coeffs

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_lucky_prime!(luckyprimes)
    @log level = -2 "The first lucky prime is $prime"
    @log level = -2 "Reducing input generators modulo $prime"

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
    @log level = -5 "Reduced coefficients are" basis_ff.coeffs

    @log level = -5 "Before F4" basis_ff
    params_zp = params_mod_p(params, prime)
    f4!(ring_ff, basis_ff, pairset, hashtable, tracer, params_zp)
    @log level = -5 "After F4:" basis_ff
    # NOTE: basis_ff may not own its coefficients, one should not mutate it
    # directly further in the code

    push!(state.gb_coeffs_ff_all, basis_ff.coeffs)

    # Reconstruct coefficients and write results to the accumulator.
    # CRT reconstrction is trivial here.
    @log level = -2 "Reconstructing coefficients from Z_$prime to QQ"
    if crt_algorithm === :incremental
        full_incremental_crt_reconstruct!(state, luckyprimes)
    else
        full_simultaneous_crt_reconstruct!(state, luckyprimes)
    end
    success_reconstruct = full_rational_reconstruct!(state, luckyprimes)
    @log level = -5 "Reconstructed coefficients" state.gb_coeffs_qq
    @log level = -2 "Successfull reconstruction: $success_reconstruct"

    correct_basis = false
    if success_reconstruct
        @log level = -2 "Verifying the correctness of reconstruction"
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
        @log level = -2 "Passed correctness check: $correct_basis"
        # At this point, if the constructed basis is correct, we return it.
        if correct_basis
            # take monomials from the basis modulo a prime
            gb_monoms, _ = export_basis_data(basis_ff, hashtable)
            # take coefficients from the reconstrcted basis
            gb_coeffs_qq = state.gb_coeffs_qq
            return gb_monoms, gb_coeffs_qq
        end
    end

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    batchsize = 1
    batchsize_multiplier = 1.4
    @log level = -2 """
      Preparing to compute bases in batches.. 
      The initial size of the batch is $batchsize. 
      The size increases in a geometric progression. 
      The batch size multiplier is $batchsize_multiplier.
      """
    if !tracer.ready_to_use
        @log level = -2 """
          The tracer is disabled until the shape of the basis is not determined via majority vote.
          The threshold for the majority vote is $(params.majority_threshold)
          """
    end

    # CRT and rational reconstrction settings
    indices_selection = Vector{Tuple{Int, Int}}(undef, length(state.gb_coeffs_zz))
    k = 1
    for i in 1:length(state.gb_coeffs_zz)
        l = length(state.gb_coeffs_zz[i])
        nl = max(isqrt(l) - 1, 1)
        if isone(l)
            continue
        end
        for j in 1:nl
            if k > length(indices_selection)
                resize!(indices_selection, 2 * length(indices_selection))
            end
            indices_selection[k] = (i, rand(2:l))
            k += 1
        end
    end
    resize!(indices_selection, k - 1)
    unique!(indices_selection)
    @log level = -2 "" length(indices_selection) length(state.gb_coeffs_zz) sum(
        map(length, state.gb_coeffs_zz)
    )

    if crt_algorithm === :incremental
        partial_incremental_crt_reconstruct!(state, luckyprimes, indices_selection)
    else
        partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)
    end

    iters = 0
    while !correct_basis
        @log level = -2 "Iteration # $iters of modular Groebner"

        for j in 1:batchsize
            prime = next_lucky_prime!(luckyprimes)
            @log level = -3 "The lucky prime is $prime"
            @log level = -3 "Reducing input generators modulo $prime"

            # Perform reduction modulo prime and store result in basis_ff
            ring_ff, basis_ff =
                reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
            params_zp = params_mod_p(params, prime)

            f4!(ring_ff, basis_ff, pairset, hashtable, tracer, params_zp)

            push!(state.gb_coeffs_ff_all, basis_ff.coeffs)

            if !majority_vote!(state, basis_ff, tracer, params)
                @log level = -3 "Majority vote is not conclusive, aborting reconstruction!"
                continue
            end

            # @log level = -2 "Reconstructing coefficients using CRT"
            # @log level = -4 "Reconstructing coefficients from Z_$(luckyprimes.modulo) * Z_$(prime) to Z_$(luckyprimes.modulo * prime)"
            if crt_algorithm === :incremental
                partial_incremental_crt_reconstruct!(state, luckyprimes, indices_selection)
            end
        end
        # @log level = -2 "Reconstructing coefficients using CRT"
        if crt_algorithm === :simultaneous
            partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)
        end

        @log level = -2 "Reconstructing coefficients to QQ"
        @log level = -4 "Reconstructing coefficients from Z_$(luckyprimes.modulo * prime) to QQ"
        success_reconstruct =
            partial_rational_reconstruct!(state, luckyprimes, indices_selection)
        @log level = -2 "Reconstruction successfull: $success_reconstruct"
        @log level = -2 """
          Used $(length(luckyprimes.primes)) primes in total over $(iters + 1) iterations.
          The current batch size is $batchsize.
          """

        # @info "" state.gb_coeffs_ff_all
        # @info "" state.selected_coeffs_zz
        # @info "" state.selected_prev_coeffs_zz
        # @info "" state.selected_coeffs_qq

        if !success_reconstruct
            @log level = -2 "Partial rational reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                @log level = -2 "Heuristic check failed for partial reconstruction"
                iters += 1
                batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
                continue
            end
        end

        # Perform full reconstruction
        @log level = -2 "Performing full CRT and rational reconstruction.."
        if crt_algorithm === :incremental
            full_incremental_crt_reconstruct!(state, luckyprimes)
        else
            full_simultaneous_crt_reconstruct!(state, luckyprimes)
        end
        success_reconstruct = full_rational_reconstruct!(state, luckyprimes)

        if !success_reconstruct
            @log level = -2 "Full reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
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
        batchsize = get_next_batchsize(batchsize, batchsize_multiplier)
    end

    @log level = -2 "Correctness check passed!"
    @log level = -2 "Used $(length(luckyprimes.primes)) primes in total over $(iters) iterations"

    # Construct the output basis.
    # Take monomials from the basis modulo a prime
    gb_monoms, _ = export_basis_data(basis_ff, hashtable)
    # Take coefficients from the reconstructed basis
    gb_coeffs_qq = state.gb_coeffs_qq

    return gb_monoms, gb_coeffs_qq
end
