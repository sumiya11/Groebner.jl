
"""
    struct Statistics


"""
mutable struct Statistics
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
