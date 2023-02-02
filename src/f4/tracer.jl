
mutable struct Tracer
    is_iteration_redundant::Vector{Int}
    pairset_size::Int
    basis_ntotal::Int
    isready::Bool

    function Tracer()
        new(Vector{Int}(undef, 0), 0, 0, false)
    end
end

isready(tr::Tracer) = tr.isready
is_iteration_redundant(tr::Tracer, i::Integer) = @inbounds tr.is_iteration_redundant[i] != 0
final_basis_size(tr::Tracer) = tr.basis_ntotal

function set_ready!(tr::Tracer)
    if !isready(tr)
        tr.isready = true
    end
    nothing
end

function set_final_basis!(tr::Tracer, basis_size::Integer)
    if !isready(tr)
        tr.basis_ntotal = basis_size
    end
    nothing
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
