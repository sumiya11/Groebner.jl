
@testset "SI comparator bug" begin
    # this will crash if the comparator is invalid
    function bug(gens::Vector{Int}, exps::Vector{Vector{T}}) where {T}
        inds = collect(1:length(gens))
        cmp  = (x, y) -> Groebner.monom_isless(exps[gens[x]],exps[gens[y]], Groebner.DegRevLex())
        sort!(inds, lt=cmp)
        @test true
    end

    gens = [2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 
            2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2]

    exps = Groebner.PowerVector{UInt64}[
        [0x0000, 0x0000], [0x0002, 0x0002], [0x0002, 0x0002]
    ]

    bug(gens, exps)
end
