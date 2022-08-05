
#=
    Same as f4.jl, but with tracing.
=#

mutable struct Tracer
    isredundant_iter::Vector{Int}
    pairset_size::Int
    basis_ntotal::Int
    ready::Bool

    function Tracer()
        new(Vector{Int}(undef, 0), 0, 0, false)
    end
end