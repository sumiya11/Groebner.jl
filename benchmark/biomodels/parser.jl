
import AbstractAlgebra

function read_BIOMDs(nspecies)
    @info "Reading biomodels"
    systems = []
    for smth in readdir(abspath("benchmark/biomodels"))
        system = []

        @info "" smth
        !(occursin("BIOMD", smth)) && continue

        io = open(abspath("benchmark/biomodels/$smth/species_map.txt"), "r")
        vs = map(strip ∘ first, map(split, readlines(io)))
        close(io)
        # @info "variables" vs
        !(length(vs) in nspecies) && continue

        io = open(abspath("benchmark/biomodels/$smth/parameters.txt"), "r")
        params = readlines(io)
        params = map(p -> replace(p, "/"=>"//"), params)
        params = map(params) do p
            p = split(p, " = ")
            if occursin("//", p[2])
                n, d = split(p[2], "//")
                "$(p[1]) = BigInt($n)//BigInt($d)"
            else
                "$(p[1]) = BigInt($(p[2]))"
            end
        end
        close(io)

        io = open(abspath("benchmark/biomodels/$smth/odes.txt"), "r")
        odes = readlines(io)
        close(io)

        symvs = map(x -> Symbol(x), vs)
        expvs = Meta.parse(join(symvs, ","))
        R, xs = eval(:((R, $expvs) = AbstractAlgebra.PolynomialRing(AbstractAlgebra.QQ, $symvs)))

        for p in params
            p = eval(Meta.parse(p))
        end

        for ode in odes
            ode = replace(ode, "{"=>"","}"=>"", "/"=>"//",","=>" ")
            ode = strip(last(split(ode, "=")))
            poly = eval(Meta.parse(ode))
            push!(system, poly)
        end

        system = map(system) do f
            if f isa AbstractAlgebra.Generic.MPoly
                f
            elseif f isa AbstractAlgebra.Generic.Frac
                lc = AbstractAlgebra.leading_coefficient(denominator(f))
                AbstractAlgebra.map_coefficients(c -> c // lc, numerator(f))
            else
                R(f)
            end
        end
        system = filter(!iszero, system)
        push!(systems, (smth, map(R, system)))
    end

    @info "Loaded" length(systems)
    systems
end
