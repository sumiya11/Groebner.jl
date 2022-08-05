
function exponent_isless_drl(ea::Vector{T}, eb::Vector{T}) where {T}
    if ea[end] < eb[end]
        return true
    elseif ea[end] != eb[end]
        return false
    end

    i = length(ea) - 1
    while i > 1 && ea[i] == eb[i]
        i -= 1
    end

    return ea[i] < eb[i] ? true : false
end

function bug(gens::Vector{Int}, exps::Vector{Vector{T}}) where {T}

    inds = collect(1:length(gens))

    cmp  = (x, y) -> (
                    println("x=$x y=$y"); 
                    exponent_isless_drl(exps[gens[x]],exps[gens[y]])
                    )

    sort!(inds, lt=cmp)
end

gens = [2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2]

exps = Vector{UInt16}[[0x0000, 0x0000], [0x0002, 0x0002], [0x0002, 0x0002]]

bug(gens, exps)
