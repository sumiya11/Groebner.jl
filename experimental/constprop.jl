module Foo

vanish_enabled() = true

macro vanish(expr)
    m = __module__
    quote
        if !$m.vanish_enabled()
            $(esc(expr))
        else
            nothing
        end
    end
end

function critical_loop(x)
    for i in eachindex(x)
        @vanish println(x[i])
        x[i] = x[i]^2
    end
end

end
