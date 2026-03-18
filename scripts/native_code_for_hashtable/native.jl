using InteractiveUtils, Groebner, Random

function get_native_code(N, MonomT)
    ord = Groebner.DegRevLex()
    ring = Groebner.PolyRing(N, ord, Groebner.CoeffModular(2), :zp)
    rng = Random.default_rng()

    ht = Groebner.hashtable_initialize(ring, rng, MonomT)

    monom_str = string(MonomT)
    monom_name = replace(monom_str, "Groebner." => "", "{" => "_", "}" => "", ", " => "_", " " => "")
    
    func_name = "hashtable_insert"
    filename = "native_$(monom_name)_N$(N)_$(func_name).txt"
    mkpath(joinpath(@__DIR__, "native"))
    filepath = joinpath(@__DIR__, "native", filename)
    
    @info "Generating native code for $monom_name with $N vars to $filepath"
    
    open(filepath, "w") do io
        code_native(
            io,
            Groebner.hashtable_insert!,
            Tuple{
                typeof(ht),
                MonomT
            },
            # debuginfo=:none
        )
    end
end

N = 11

configs = [
    (N, Groebner.PackedTuple2{UInt64, UInt8}),
    (N, Vector{UInt8}),
    (N, Groebner.FixedVector{nextpow(2, N), UInt8}),
    (N, Groebner.FixedMonom{nextpow(2, N), UInt8}),
    (N, Groebner.FixedMonomNoDeg{nextpow(2, N), UInt8}),
    (N, Groebner.NibbleMonom{nextpow(2, N) ÷ 2})
]

for (N, MonomT) in configs
    get_native_code(N, MonomT)
end
