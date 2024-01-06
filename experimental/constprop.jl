module Foo

macro unreachable()
    :(@assert false)
end

function do_stuff(x)
    if x > 0
        x^2
    else
        @unreachable
        x
    end
end

end
