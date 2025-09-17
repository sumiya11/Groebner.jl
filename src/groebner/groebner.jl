# This file is a part of Groebner.jl. License is GNU GPL v2.

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
    ring, monoms, coeffs = ir_ensure_valid(ring, monoms, coeffs)
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
            @info """
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
    _, ring2, monoms2, coeffs2 = ir_convert_ir_to_internal(ring, monoms, coeffs, params)
    gb_monoms2, gb_coeffs2 = groebner2(ring2, monoms2, coeffs2, params)
    gb_monoms, gb_coeffs = ir_convert_internal_to_ir(ring2, gb_monoms2, gb_coeffs2, params)
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
        _, ring, monoms, coeffs, params = homogenize_generators(ring, monoms, coeffs, params)
    end

    gb_monoms, gb_coeffs = _groebner2(ring, monoms, coeffs, params)

    if params.homogenize
        ring, gb_monoms, gb_coeffs, params =
            dehomogenize_generators(ring, gb_monoms, gb_coeffs, params)
        if params.reduced
            gb_monoms, gb_coeffs = _groebner2(ring, gb_monoms, gb_coeffs, params)
        end
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
    basis, pairset, hashtable = f4_initialize_structs(ring, monoms, coeffs, params)
    f4!(ring, basis, pairset, hashtable, params)
    gbmonoms, gbcoeffs = basis_export_data(basis, hashtable)
    gbmonoms, gbcoeffs
end

###
# Groebner basis over generic field. Calls F4 directly.

function _groebner2(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffGeneric}
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
        _groebner_learn_and_apply(ring, monoms, coeffs, params)
    else
        @assert params.modular_strategy === :classic_modular
        _groebner_classic_modular(ring, monoms, coeffs, params)
    end
end

###
# Learn & Apply startegy

# The next batch is a multiple of the previous ones
function get_next_batch(
    primes_used::Int,
    prev_batch::Int,
    batch_scaling::Float64,
    composite::Int,
    tasks::Int
)
    new_batch = max(composite, round(Int, primes_used * batch_scaling))
    new_batch = align_up(new_batch, composite * tasks)
    max(new_batch, prev_batch)
end

@timeit _TIMER function _groebner_guess_lucky_prime(
    state::ModularState,
    ring::PolyRing,
    basis_zz::Basis,
    pairset::Pairset,
    hashtable::MonomialHashtable{M},
    params::AlgorithmParameters
) where {M <: Monom}
    prime_1 = modular_random_prime(state, params.rng)
    ring_ff_1, basis_ff_1 = modular_reduce_mod_p!(ring, basis_zz, prime_1, deepcopy=true)
    params_zp = param_mod_p(params, prime_1)
    f4!(ring_ff_1, basis_ff_1, pairset, hashtable, params_zp)

    prime_2 = modular_random_prime(state, params.rng)
    ring_ff_2, basis_ff_2 = modular_reduce_mod_p!(ring, basis_zz, prime_2, deepcopy=true)
    params_zp = param_mod_p(params, prime_2)
    f4!(ring_ff_2, basis_ff_2, pairset, hashtable, params_zp)

    if basis_ff_1.monoms == basis_ff_2.monoms
        return prime_1
    end

    prime_3 = modular_random_prime(state, params.rng)
    ring_ff_3, basis_ff_3 = modular_reduce_mod_p!(ring, basis_zz, prime_3, deepcopy=true)
    params_zp = param_mod_p(params, prime_3)
    f4!(ring_ff_3, basis_ff_3, pairset, hashtable, params_zp)

    if basis_ff_1.monoms == basis_ff_3.monoms
        return prime_1
    elseif basis_ff_2.monoms == basis_ff_3.monoms
        return prime_2
    else
        throw("Cannot find a lucky prime. Groebner.jl gives up.")
    end

    prime_1
end

###
# Learn & Apply strategy

# internal function for multi-threading
function _process_chunk(
    chunk::Vector{Int},
    ring::PolyRing,
    basis_zz::Basis,
    composite::Int,
    trace::Trace,
    primes::Vector{Int32},
    params::AlgorithmParameters,
    gb_coeffs::Vector{Tuple{Int32, Vector{Vector{Int32}}}},
)
    for i in 1:length(chunk)
        primes_x = ntuple(k -> primes[(i - 1) * composite + k], composite)
        ring_zp_x, basis_zp_x = modular_reduce_mod_p_in_batch!(ring, basis_zz, primes_x)
        params_zp = param_mod_p(
            params,               # can be mutated later
            CompositeNumber{composite, Int32}(primes_x),
            using_wide_type_for_coeffs=false
        )
        trace.buf_basis = basis_zp_x
        trace.ring = ring_zp_x
        flag = f4_apply!(trace, params_zp)
        !flag && return
        gb_coeffs_unpacked = ir_unpack_composite_coefficients(trace.gb_basis.coeffs)
        append!(gb_coeffs, collect(zip(primes_x, gb_coeffs_unpacked)))
    end
    return
end

@timeit _TIMER2 function _groebner_learn_and_apply(
    ring::PolyRing,
    monoms::Vector{Vector{M}},
    coeffs::Vector{Vector{C}},
    params::AlgorithmParameters
) where {M <: Monom, C <: CoeffQQ}
    # Guess lucky prime
    basis, pairset, hashtable, permutation =
        f4_initialize_structs(ring, monoms, coeffs, params, make_monic=false)
    basis_zz = clear_denominators!(basis, deepcopy=false)
    state = ModularState{BigInt, C, Int32}(basis_zz.coeffs)
    @timeit _TIMER2 "_groebner_guess_lucky_prime" prime = _groebner_guess_lucky_prime(state, ring, basis_zz, pairset, hashtable, params)

    # Learn
    ring_ff, basis_ff = modular_reduce_mod_p!(ring, basis_zz, prime, deepcopy=true)
    params_zp = param_mod_p(params, prime)
    trace = trace_initialize(ring_ff, basis_ff, hashtable, permutation, params_zp)
    @timeit _TIMER2 "f4_learn!" f4_learn!(trace, pairset, params_zp)

    # CRT and rational reconstruction settings
    # TODO: no need to deepcopy!
    push!(state.used_primes, prime)
    push!(state.gb_coeffs_ff_all, deepcopy(trace.gb_basis.coeffs))
    modular_prepare!(state)
    witness_set = modular_witness_set(state.gb_coeffs_zz, params)

    # Continue to compute Groebner bases modulo different primes in batches. 
    composite = params.composite
    tasks = params.tasks
    batch_scaling = params.batch_scaling
    primes_used = 1
    correct_basis = false
    batch = composite * tasks

    # Initialize a tracer that computes the bases using composite numbers
    trace_x = trace_copy(trace, PolynomialRepresentation(M, CompositeNumber{composite, Int32}, false))

    # Task local buffers
    task_results = Vector{Task}(undef, tasks)
    tasklocal_trace_x = map(_ -> deepcopy(trace_x), 1:tasks)
    tasklocal_gb_coeffs = Vector{Vector{Tuple{Int32, Vector{Vector{Int32}}}}}(undef, tasks)
    for i in 1:tasks
        tasklocal_gb_coeffs[i] = Vector{Tuple{Int, Vector{Vector{Int32}}}}()
    end
    tasklocal_params = map(_ -> deepcopy(params), 1:tasks)

    while !correct_basis
        batch = get_next_batch(primes_used, batch, batch_scaling, composite, tasks)
        @invariant iszero(batch % composite) && iszero(batch % tasks)

        tasklocal_primes =
            [map(_ -> Int32(modular_next_prime!(state)), 1:div(batch, tasks)) for _ in 1:tasks]
        for i in 1:tasks
            empty!(tasklocal_gb_coeffs[i])
        end

        data_chunks = split_round_robin(1:composite:batch, tasks)
        for (tid, chunk) in enumerate(data_chunks)
            task = @spawn begin
                _process_chunk(
                    chunk,
                    ring,
                    basis_zz,
                    composite,
                    tasklocal_trace_x[tid],
                    tasklocal_primes[tid],
                    tasklocal_params[tid],
                    tasklocal_gb_coeffs[tid]
                )
            end
            task_results[tid] = task
        end
        @timeit _TIMER2 "f4_apply!" for task in task_results
            wait(task)
        end

        primes_used += batch

        tasklocal_gb_coeffs_union = reduce(vcat, tasklocal_gb_coeffs)
        sort!(tasklocal_gb_coeffs_union, by=first, rev=true)
        for (prime_, coeffs_ff_) in tasklocal_gb_coeffs_union
            push!(state.used_primes, prime_)
            push!(state.gb_coeffs_ff_all, coeffs_ff_)
        end

        @timeit _TIMER2 "crt_vec_partial!" crt_vec_partial!(
            state.gb_coeffs_zz,
            state.modulo,
            state.gb_coeffs_ff_all,
            state.used_primes,
            witness_set,
            state.crt_mask
        )
        @timeit _TIMER2 "ratrec_vec_partial!" success_reconstruct = ratrec_vec_partial!(
            state.gb_coeffs_qq,
            state.gb_coeffs_zz,
            state.modulo,
            witness_set,
            state.ratrec_mask
        )

        if !success_reconstruct
            continue
        end

        if params.heuristic_check
            success_check =
                modular_lift_heuristic_check_partial(state.gb_coeffs_qq, state.modulo, witness_set)
            if !success_check
                continue
            end
        end

        # Perform full reconstruction
        @timeit _TIMER2 "crt_vec_full!" crt_vec_full!(
            state.gb_coeffs_zz,
            state.modulo,
            state.gb_coeffs_ff_all,
            state.used_primes,
            state.crt_mask;
            tasks=params.tasks
        )
        @timeit _TIMER2 "ratrec_vec_full!" success_reconstruct = ratrec_vec_full!(
            state.gb_coeffs_qq,
            state.gb_coeffs_zz,
            state.modulo,
            state.ratrec_mask,
            tasks=params.tasks
        )

        # This should happen rarely
        if !success_reconstruct
            continue
        end

        @timeit _TIMER2 "modular_lift_check!" correct_basis =
            modular_lift_check!(state, ring_ff, basis, basis_zz, trace.gb_basis, hashtable, params)
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
    basis, pairset, hashtable =
        f4_initialize_structs(ring, monoms, coeffs, params, make_monic=false)
    basis_zz = clear_denominators!(basis, deepcopy=false)
    state = ModularState{BigInt, C, CoeffModular}(basis_zz.coeffs)
    prime = _groebner_guess_lucky_prime(state, ring, basis_zz, pairset, hashtable, params)

    ring_ff, basis_ff = modular_reduce_mod_p!(ring, basis_zz, prime, deepcopy=true)
    params_zp = param_mod_p(params, prime)
    f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)
    # NOTE: basis_ff may not own its coefficients, one should not mutate it
    # directly further in the code

    # CRT and rational reconstruction settings
    push!(state.used_primes, prime)
    push!(state.gb_coeffs_ff_all, basis_ff.coeffs)
    modular_prepare!(state)
    witness_set = modular_witness_set(state.gb_coeffs_zz, params)

    # At this point, either the reconstruction or the correctness check failed.
    # Continue to compute Groebner bases modulo different primes in batches. 
    composite = 1
    tasks = 1
    batch_scaling = params.batch_scaling
    primes_used = 1
    correct_basis = false
    batch = 1

    while !correct_basis
        batch = get_next_batch(primes_used, batch, batch_scaling, composite, tasks)

        for j in 1:batch
            prime = modular_next_prime!(state)
            ring_ff, basis_ff = modular_reduce_mod_p!(ring, basis_zz, prime, deepcopy=true)
            params_zp = param_mod_p(params, prime)
            f4!(ring_ff, basis_ff, pairset, hashtable, params_zp)
            push!(state.used_primes, prime)
            push!(state.gb_coeffs_ff_all, basis_ff.coeffs)
            primes_used += 1
        end

        crt_vec_partial!(
            state.gb_coeffs_zz,
            state.modulo,
            state.gb_coeffs_ff_all,
            state.used_primes,
            witness_set,
            state.crt_mask
        )
        success_reconstruct = ratrec_vec_partial!(
            state.gb_coeffs_qq,
            state.gb_coeffs_zz,
            state.modulo,
            witness_set,
            state.ratrec_mask
        )

        if !success_reconstruct
            continue
        end

        if params.heuristic_check
            success_check =
                modular_lift_heuristic_check_partial(state.gb_coeffs_qq, state.modulo, witness_set)
            if !success_check
                continue
            end
        end

        # Perform full reconstruction
        crt_vec_full!(
            state.gb_coeffs_zz,
            state.modulo,
            state.gb_coeffs_ff_all,
            state.used_primes,
            state.crt_mask,
            tasks=params.tasks
        )
        success_reconstruct = ratrec_vec_full!(
            state.gb_coeffs_qq,
            state.gb_coeffs_zz,
            state.modulo,
            state.ratrec_mask,
            tasks=params.tasks
        )

        if !success_reconstruct
            continue
        end

        correct_basis =
            modular_lift_check!(state, ring_ff, basis, basis_zz, basis_ff, hashtable, params)
    end

    gb_monoms, _ = basis_export_data(basis_ff, hashtable)
    gb_coeffs_qq = state.gb_coeffs_qq

    return gb_monoms, gb_coeffs_qq
end
