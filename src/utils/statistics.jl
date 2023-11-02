# This file is WIP and is not used at the moment

mutable struct Statistics
    # Core stats -- F4, matrix reduction
    a::Any

    # Modular computation
    # (is not printed when computing over finite fields)
    b::Any

    # Is printed sometimes
    misc::Any

    info::Dict{Symbol, Any}
end

macro record(kws)
    @assert kws.head == :(=)
    key, value = kws.args[1], kws.args[2]
    file, line = String(__source__.file), Int(__source__.line)
    esc(_record(file, line, key, value))
end

const _global_statistics = Statistics(1, 1, 2, Dict{Symbol, Any}())
# const _perf_id_to_index = Dict{Any, Int}()

function _record(file, line, key, value)
    key_ = Meta.quot(key)
    body = quote
        result = $value
        $(@__MODULE__)._global_statistics.info[$key_] = result
        result
    end
    quote
        if $(@__MODULE__).logging_enabled()
            $body
        else
            nothing
        end
    end
end
