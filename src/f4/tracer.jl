
# Tracer
# Contains some information about previous finite field runs of F4
# within the current modular run, which could be useful for subsequent rungs.
mutable struct Tracer
    is_iteration_redundant::Vector{Int}
    pairset_size::Int
    basis_ntotal::Int
    ready_to_use::Bool
    function Tracer()
        new(Vector{Int}(undef, 0), 0, 0, false)
    end
    function Tracer(params)
        Tracer()
    end
end

isready(tr::Tracer) = tr.ready_to_use
is_iteration_redundant(tr::Tracer, i::Integer) = @inbounds tr.is_iteration_redundant[i] != 0
final_basis_size(tr::Tracer) = tr.basis_ntotal

function set_ready!(tr::Tracer)
    tr.ready_to_use = true
end

function set_final_basis!(tr::Tracer, basis_size::Integer)
    tr.basis_ntotal = basis_size
end

function update_tracer_pairset!(tr::Tracer, pairset_size::Integer)
    if !isready(tr)
        tr.pairset_size = max(pairset_size, tr.pairset_size)
    end
    nothing
end

function update_tracer_iteration!(tr::Tracer, isredundant::Bool)
    if !isready(tr)
        push!(tr.is_iteration_redundant, 0)
        if isredundant
            tr.is_iteration_redundant[end] = 1
        end
    end
    nothing
end
