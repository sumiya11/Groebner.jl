# Very lightweight tracing in f4 (for the real stuff, see trace.jl)

mutable struct TinyTraceF4
    is_iteration_redundant::Vector{Int}
    pairset_size::Int
    basis_nfilled::Int
    ready_to_use::Bool

    TinyTraceF4(params) = TinyTraceF4()

    function TinyTraceF4()
        new(Vector{Int}(), 0, 0, false)
    end
end

isready(tr::TinyTraceF4) = tr.ready_to_use
is_iteration_redundant(tr::TinyTraceF4, i::Integer) =
    @inbounds tr.is_iteration_redundant[i] != 0
final_basis_size(tr::TinyTraceF4) = tr.basis_nfilled

function set_ready!(tr::TinyTraceF4)
    tr.ready_to_use = true
end

function set_final_basis!(tr::TinyTraceF4, basis_size::Integer)
    tr.basis_nfilled = basis_size
end

function update_tracer_pairset!(tr::TinyTraceF4, pairset_size::Integer)
    if !isready(tr)
        tr.pairset_size = max(pairset_size, tr.pairset_size)
    end
    nothing
end

function update_tracer_iteration!(tr::TinyTraceF4, isredundant::Bool)
    if !isready(tr)
        push!(tr.is_iteration_redundant, 0)
        if isredundant
            tr.is_iteration_redundant[end] = 1
        end
    end
    nothing
end
