
function f_continue(x)
    i = 1
    @inbounds while i < length(x)
        if iszero(x[i])
            i += 1
            continue
        end
        i += 2
    end
    i
end

function f_goto(x)
    i = 1
    @label Letsgo
    @inbounds while i < length(x)
        if iszero(x[i])
            i += 1
            @goto Letsgo
        end
        i += 2
    end
    i
end

using BenchmarkTools

x = rand(Int, 1000)
x[rand(1:1000, 500)] .= 0

@benchmark f_continue($x)

@benchmark f_goto($x)
