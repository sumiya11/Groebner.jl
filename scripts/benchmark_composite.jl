using Pkg;
Pkg.develop(path=joinpath((@__DIR__), "../"))
# Pkg.add(url="https://github.com/sumiya11/Groebner.jl", rev="benchmark-composite")

using Plots, Groebner, BenchmarkTools

SAVE_DIR = (@__DIR__)
BOOT = 10
COMPOSITE = (1, 2, 4, 8, 16, 32, 64)

@info "Benchmarking composite widths: $COMPOSITE"
@info "Using BOOT=$BOOT"
@info "Saving result to $SAVE_DIR"

systems = [
    ("chandra-9", Groebner.Examples.chandran(9)),
    ("chandra-10", Groebner.Examples.chandran(10)),
    ("cyclic-7", Groebner.Examples.cyclicn(7)),
    ("cyclic-8", Groebner.Examples.cyclicn(8)),
    ("eco-11", Groebner.Examples.eco11()),
    ("katsura-9", Groebner.Examples.katsuran(9)),
    ("katsura-10", Groebner.Examples.katsuran(10)),
    ("hexapod", Groebner.Examples.hexapod()),
    ("hiv2", Groebner.Examples.HIV2()),
    ("alea6", Groebner.Examples.alea6())
]

data = Dict()
for (name, sys) in systems
    @info "Running $name.."
    data[name] = Dict("time" => [], "memory" => [], "primes" => 0)
    groebner(sys, threaded=:no)
    data[name]["primes"] = Groebner.PRIMES[]
    for _composite in COMPOSITE
        ts = [@elapsed groebner(sys; _composite=_composite, threaded=:no) for _ in 1:BOOT]
        a = @allocated groebner(sys; _composite=_composite, threaded=:no)
        t = minimum(ts)
        push!(data[name]["time"], t)
        push!(data[name]["memory"], a)
    end
    data[name]["time_ratio"] = data[name]["time"] / data[name]["time"][1]
    data[name]["memory_ratio"] = data[name]["memory"] / data[name]["memory"][1]
end

function draw()
    begin
        dpi = 200
        figsize = (700, 900)
        linewidth = 2
        titlesize = 8
        linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
        p1 = plot(
            size=figsize,
            xticks=(1:length(COMPOSITE), collect(COMPOSITE)),
            yticks=[0.25, 0.5, 0.75, 1.0, 1.25, 1.5],
            ylims=(0.30, 1.70),
            xlabel="composite width",
            ylabel="time ratio",
            legend=:topleft,
            title="Runtime",
            titlesize=titlesize
        )
        p2 = plot(
            size=figsize,
            xticks=(1:length(COMPOSITE), collect(COMPOSITE)),
            yticks=[0.25, 0.5, 0.75, 1.0, 1.25, 1.5],
            ylims=(0.20, 1.60),
            xlabel="composite width",
            ylabel="memory ratio",
            legend=:topleft,
            title="Memory",
            titlesize=titlesize
        )
        for (i, (name, _)) in enumerate(systems)
            name_padded = rpad(name, 13)
            label1 = "$name_padded ($(round(data[name]["time"][1], digits=1)) s)"
            plot!(
                p1,
                1:length(COMPOSITE),
                data[name]["time_ratio"],
                label=label1,
                alpha=0.7,
                linewidth=linewidth,
                linestyle=linestyles[i % length(linestyles) + 1]
            )
            primes = data[name]["primes"]
            annotate!(
                p1,
                3.0,
                1.6 - 0.048 * (i - 1),
                Plots.text(
                    "∘ $(name) - $(primes) primes",
                    halign=:left,
                    color="black",
                    pointsize=9
                )
            )

            label2 = "$name_padded ($(BenchmarkTools.prettymemory(data[name]["memory"][1])))"
            plot!(
                p2,
                1:length(COMPOSITE),
                data[name]["memory_ratio"],
                label=label2,
                alpha=0.7,
                linewidth=linewidth,
                linestyle=linestyles[i % length(linestyles) + 1]
            )
        end
        P = plot(p1, p2, layout=(2, 1))
        savefig(P, joinpath([SAVE_DIR, "plot.pdf"]))
        P
    end
end

draw()
