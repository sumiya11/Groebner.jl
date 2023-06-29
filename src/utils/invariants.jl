
"""
    macro invariant

1
"""
macro invariant(arg)
    dir = @__DIR__
    file, line = String(__source__.file), Int(__source__.line)
    esc(:(
        if $(@__MODULE__).invariants_enabled()
            _invariant($dir, $file, $line, $arg, $(string(arg)))
        else
            nothing
        end
    ))
end

function _invariant(dir, file, line, arg, str)
    if !arg
        throw(AssertionError("$str"))
    end
    nothing
end
