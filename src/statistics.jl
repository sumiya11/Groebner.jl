
"""
    struct Statistics


"""
mutable struct Statistics
    # Runtimes spent in different components of F4
    time_select::Any
    time_symbol::Any
    time_linalg::Any
    time_update::Any
    time_reduce::Any
    time_learn_trace::Any
    time_apply_trace::Any

    # Runtimes spent in high level groebner basis
    time_input_output::Any
    time_modular_reduce::Any
    time_crt_reconstruct::Any
    time_modular_reconstruct::Any
    time_correctness_check::Any

    num_iters_f4::Any
    num_pairs::Any
    num_monoms::Any
    num_generators::Any
    num_div_checks::Any
    num_divmask_hits::Any

    used_primes::Any

    # Estimated sizes in bytes
    size_matrix::Any
    size_pairset::Any
    size_basis::Any
end
