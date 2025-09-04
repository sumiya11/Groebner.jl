using Nemo, Groebner, BenchmarkTools, PrettyTables, Printf

k = Nemo.GF(2^30+3)
systems = [
    ("chandra-9", Groebner.Examples.chandran(9, k=k)),
    ("chandra-10", Groebner.Examples.chandran(10, k=k)),
    ("cyclic-7", Groebner.Examples.cyclicn(7, k=k)),
    ("eco-11", Groebner.Examples.econ(11, k=k)),
    ("noon-7", Groebner.Examples.noonn(7, k=k)),
    ("katsura-9", Groebner.Examples.katsuran(9, k=k)),
    ("katsura-10", Groebner.Examples.katsuran(10, k=k)),
    ("hexapod", Groebner.Examples.hexapod(k=k)),
]

println("Running the following systems: ", map(first, systems))

data = []
for (name, sys) in systems
    @info "Running $name.."
    t1, t2, t3, t4 = 0, 0, 0, 0
    for _ in 1:2
        t1 = @elapsed groebner(sys; threaded=:no)
        t2 = @elapsed trace, _ = groebner_learn(sys)
        t3 = @elapsed groebner_apply!(trace, sys)
        t3 = @elapsed groebner_apply!(trace, sys)
        t4 = @elapsed groebner_apply!(trace, (sys, sys, sys, sys))
        t4 = @elapsed groebner_apply!(trace, (sys, sys, sys, sys))
    end
    push!(data, [name, t1, t2, t3, t4])
end

matrix = permutedims(reduce(hcat, data))
pretty_table(
    matrix, 
    column_labels=[
        [EmptyCells(2), MultiColumn(3, "Learn & Apply")], 
        ["Name", "Monte-Carlo", "Learn", "Apply", "Apply 4x"]
    ],
    title="Timings in seconds",
    formatters=[(v,i,j) -> v isa Number ? @sprintf("%.2f", v) : v],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded);  
)
