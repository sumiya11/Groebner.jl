# This file is a part of Groebner.jl. License is GNU GPL v2.

### 
# Backend for `groebner_with_change_matrix`

function groebner_with_change_matrix0(polynomials, options)
    isempty(polynomials) && throw(DomainError("Empty input."))
    ring, monoms, coeffs, options = io_convert_polynomials_to_ir(polynomials, options)
    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs =
        _groebner_with_change_matrix1(ring, monoms, coeffs, options)
    result = io_convert_ir_to_polynomials(ring, polynomials, gb_monoms, gb_coeffs, options)
    rows = [
        io_convert_ir_to_polynomials(ring, polynomials, m, c, options) for
        (m, c) in zip(matrix_monoms, matrix_coeffs)
    ]
    matrix = collect(transpose(reduce(hcat, rows)))
    result, matrix
end

function groebner_with_change_matrix1(ring, monoms, coeffs, options)
    ring, monoms, coeffs = ir_ensure_assumptions(ring, monoms, coeffs)
    _groebner_with_change_matrix1(ring, monoms, coeffs, options)
end

function _groebner_with_change_matrix1(ring::PolyRing, monoms, coeffs, options)
    try
        params = AlgorithmParameters(ring, options)
        return __groebner_with_change_matrix1(ring, monoms, coeffs, params)
    catch err
        if isa(err, MonomialDegreeOverflow)
            @log :info """
            Possible overflow of exponent vector detected. 
            Restarting with at least 32 bits per exponent."""
            params = AlgorithmParameters(ring, options; hint=:large_exponents)
            return __groebner_with_change_matrix1(ring, monoms, coeffs, params)
        else
            # Something bad happened.
            rethrow(err)
        end
    end
end

function __groebner_with_change_matrix1(ring, monoms, coeffs, params)
    @invariant ir_is_valid(ring, monoms, coeffs)
    _, ring2, monoms2, coeffs2 =
        io_convert_ir_to_internal(ring, monoms, coeffs, params, params.representation)
    gb_monoms2, gb_coeffs2, matrix_monoms2, matrix_coeffs2 =
        groebner_with_change_matrix2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = io_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)
    matrix_monoms_and_coeffs = [
        io_convert_internal_to_ir(ring2, _matrix_monoms2, _matrix_coeffs2, params) for
        (_matrix_monoms2, _matrix_coeffs2) in zip(matrix_monoms2, matrix_coeffs2)
    ]
    matrix_monoms, matrix_coeffs =
        map(first, matrix_monoms_and_coeffs), map(last, matrix_monoms_and_coeffs)
    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs
end

function groebner_with_change_matrix2(ring, monoms, coeffs, params)
    nonzero_indices = findall(!io_iszero_coeffs, coeffs)
    zero_indices = setdiff(collect(1:length(monoms)), nonzero_indices)

    __monoms = monoms
    __coeffs = coeffs
    _monoms = filter(!isempty, monoms)
    _coeffs = filter(!isempty, coeffs)
    if isempty(_monoms)
        return [monoms[1]],
        [coeffs[1]],
        [[monoms[1] for _ in 1:length(monoms)]],
        [[coeffs[1] for _ in 1:length(coeffs)]]
    end
    monoms, coeffs = _monoms, _coeffs

    if params.homogenize
        _, ring, monoms, coeffs = homogenize_generators!(ring, monoms, coeffs, params)
    end

    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs =
        _groebner_with_change_matrix2(ring, monoms, coeffs, params)

    if params.homogenize
        ring, gb_monoms, gb_coeffs =
            dehomogenize_generators!(ring, gb_monoms, gb_coeffs, params)
    end

    if !isempty(zero_indices)
        for i in 1:length(matrix_monoms)
            _matrix_monoms_i = resize!(similar(matrix_monoms[i]), length(__monoms))
            _matrix_coeffs_i = resize!(similar(matrix_coeffs[i]), length(__monoms))
            zero_indices = setdiff(collect(1:length(__monoms)), nonzero_indices)
            _matrix_monoms_i .= [__monoms[zero_indices[1]]]
            _matrix_coeffs_i .= [__coeffs[zero_indices[1]]]
            for j in 1:length(nonzero_indices)
                _matrix_monoms_i[nonzero_indices[j]] = matrix_monoms[i][j]
                _matrix_coeffs_i[nonzero_indices[j]] = matrix_coeffs[i][j]
            end
            matrix_monoms[i] = _matrix_monoms_i
            matrix_coeffs[i] = _matrix_coeffs_i
        end
    end

    gb_monoms, gb_coeffs, matrix_monoms, matrix_coeffs
end

###
# Groebner basis over Z_p.
# Just calls f4 directly.

function _groebner_with_change_matrix2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffZp}
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    f4!(ring, basis, pairset, hashtable, params)
    gbmonoms, gbcoeffs = basis_export_data(basis, hashtable)
    matrix_monoms, matrix_coeffs =
        basis_changematrix_export(basis, hashtable, length(monoms))
    gbmonoms, gbcoeffs, matrix_monoms, matrix_coeffs
end

###
# Groebner basis over Q.

# GB over the rationals uses modular computation.
function _groebner_with_change_matrix2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    _groebner_with_change_classic_modular(ring, monoms, coeffs, params)
end

###
# Classic multi-modular strategy

function _groebner_with_change_classic_modular(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    # NOTE: we can mutate ring, monoms, and coeffs here.
    @log :misc "Backend: classic multi-modular F4 (with change matrix)"

    # Initialize supporting structs
    state = GroebnerState{BigInt, C, CoeffModular}(params)
    basis, pairset, hashtable =
        f4_initialize_structs(ring, monoms, coeffs, params, make_monic=false)

    # Scale the input coefficients to integers to speed up the subsequent search
    # for lucky primes
    @log :all "Input polynomials" basis
    @log :misc "Clearing the denominators of the input polynomials"
    basis_zz = clear_denominators!(state.buffer, basis, deepcopy=false)
    @log :all "Integer coefficients are" basis_zz.coeffs

    # Handler for lucky primes
    luckyprimes = LuckyPrimes(basis_zz.coeffs)
    prime = next_lucky_prime!(luckyprimes)
    @log :misc "The first lucky prime is $prime"
    @log :misc "Reducing input generators modulo $prime"

    # Perform reduction modulo prime and store result in basis_ff
    ring_ff, basis_ff = reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
    @log :all "Reduced coefficients are" basis_ff.coeffs

    @log :all "Before F4" basis_ff
    params_zp = params_mod_p(params, prime)
    f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)
    @log :all "After F4:" basis_ff
    # NOTE: basis_ff may not own its coefficients, one should not mutate it
    # directly further in the code

    push!(state.gb_coeffs_ff_all, basis_ff.coeffs)

    changematrix_monoms, changematrix_coeffs =
        basis_changematrix_export(basis_ff, hashtable, length(monoms))
    push!(state.changematrix_coeffs_ff_all, changematrix_coeffs)

    # Reconstruct coefficients and write results to the accumulator.
    # CRT reconstrction is trivial here.
    @log :misc "Reconstructing coefficients from Z_$prime to QQ"
    full_simultaneous_crt_reconstruct!(state, luckyprimes)
    full_simultaneous_crt_reconstruct_changematrix!(state, luckyprimes)

    success_reconstruct = full_rational_reconstruct!(state, luckyprimes, params.use_flint)
    @log :all "Reconstructed coefficients" state.gb_coeffs_qq
    @log :misc "Successfull reconstruction: $success_reconstruct"

    changematrix_success_reconstruct =
        full_rational_reconstruct_changematrix!(state, luckyprimes, params.use_flint)
    @log :misc "Successfull change matrix reconstruction: $success_reconstruct"

    correct_basis = false
    if success_reconstruct && changematrix_success_reconstruct
        @log :misc "Verifying the correctness of reconstruction"
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
        @log :misc "Passed correctness check: $correct_basis"
        # At this point, if the constructed basis is correct, we return it.
        if correct_basis
            # take monomials from the basis modulo a prime
            gb_monoms, _ = basis_export_data(basis_ff, hashtable)
            # take coefficients from the reconstrcted basis
            gb_coeffs_qq = state.gb_coeffs_qq
            changematrix_coeffs_qq = state.changematrix_coeffs_qq
            return gb_monoms, gb_coeffs_qq, changematrix_monoms, changematrix_coeffs_qq
        end
    end

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    primes_used = 1
    batchsize = 1
    batchsize_scaling = 0.10
    @log :misc """
    Preparing to compute bases in batches.. 
    The initial size of the batch is $batchsize. 
    The batch scale factor is $batchsize_scaling."""

    tot, indices_selection = gb_modular_select_indices0(state.gb_coeffs_zz, params)
    @log :misc "Partial reconstruction: $(length(indices_selection)) indices out of $tot"

    partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

    iters = 0
    while !correct_basis
        @log :misc "Iteration # $iters of modular Groebner"

        for j in 1:batchsize
            prime = next_lucky_prime!(luckyprimes)
            @log :debug "The lucky prime is $prime"
            @log :debug "Reducing input generators modulo $prime"

            # Perform reduction modulo prime and store result in basis_ff
            ring_ff, basis_ff =
                reduce_modulo_p!(state.buffer, ring, basis_zz, prime, deepcopy=true)
            params_zp = params_mod_p(params, prime)

            f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)

            push!(state.gb_coeffs_ff_all, basis_ff.coeffs)
            _, changematrix_coeffs =
                basis_changematrix_export(basis_ff, hashtable, length(monoms))
            push!(state.changematrix_coeffs_ff_all, changematrix_coeffs)

            if !majority_vote!(state, basis_ff, params)
                @log :debug "Majority vote is not conclusive, aborting reconstruction!"
                continue
            end

            primes_used += 1
        end
        partial_simultaneous_crt_reconstruct!(state, luckyprimes, indices_selection)

        @log :misc "Partially reconstructing coefficients to QQ"
        @log :debug "Partially reconstructing coefficients from Z_$(luckyprimes.modulo * prime) to QQ"
        success_reconstruct = partial_rational_reconstruct!(
            state,
            luckyprimes,
            indices_selection,
            params.use_flint
        )
        @log :misc "Partial reconstruction successfull: $success_reconstruct"
        @log :misc """
          Used $(length(luckyprimes.primes)) primes in total over $(iters + 1) iterations.
          The current batch size is $batchsize.
          """

        if !success_reconstruct
            @log :misc "Partial rational reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        if params.heuristic_check
            success_check =
                heuristic_correctness_check(state.selected_coeffs_qq, luckyprimes.modulo)
            if !success_check
                @log :misc "Heuristic check failed for partial reconstruction"
                iters += 1
                batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
                continue
            end
        end

        # Perform full reconstruction
        @log :misc "Performing full CRT and rational reconstruction.."
        full_simultaneous_crt_reconstruct!(state, luckyprimes)

        success_reconstruct =
            full_rational_reconstruct!(state, luckyprimes, params.use_flint)

        if !success_reconstruct
            @log :misc "Full reconstruction failed"
            iters += 1
            batchsize = get_next_batchsize(primes_used, batchsize, batchsize_scaling)
            continue
        end

        full_simultaneous_crt_reconstruct_changematrix!(state, luckyprimes)
        success_reconstruct_changematrix =
            full_rational_reconstruct_changematrix!(state, luckyprimes, params.use_flint)
        if !success_reconstruct_changematrix
            @log :misc "Failed to reconstruct the change matrix"
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

    @log :misc "Correctness check passed!"
    @log :misc "Used $(length(luckyprimes.primes)) primes in total over $(iters) iterations"

    # Construct the output basis.
    # Take monomials from the basis modulo a prime
    gb_monoms, _ = basis_export_data(basis_ff, hashtable)
    # Take coefficients from the reconstructed basis
    gb_coeffs_qq = state.gb_coeffs_qq
    changematrix_coeffs_qq = state.changematrix_coeffs_qq

    return gb_monoms, gb_coeffs_qq, changematrix_monoms, changematrix_coeffs_qq
end
