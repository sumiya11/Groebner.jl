
@testset "SI comparator bug" begin

    function bug(gens::Vector{Int}, exps::Vector{Vector{T}}) where {T}

        inds = collect(1:length(gens))

        cmp  = (x, y) -> Groebner.exponent_isless_drl(exps[gens[x]],exps[gens[y]])

        sort!(inds, lt=cmp)
    end

    gens = [2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 
            2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2]

    exps = Vector{UInt16}[[0x0000, 0x0000], [0x0002, 0x0002], [0x0002, 0x0002]]

    bug(gens, exps)

end