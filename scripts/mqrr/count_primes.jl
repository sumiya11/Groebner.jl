using AbstractAlgebra, Groebner, BenchmarkTools, PrettyTables, Printf, TimerOutputs

include("gizmos.jl")

k = AbstractAlgebra.GF(2^30+3)
system_solving = [
    ("chandra-9", Groebner.Examples.chandran(9, k=k)),
    ("chandra-10", Groebner.Examples.chandran(10, k=k)),
    ("chandra-11", Groebner.Examples.chandran(11, k=k)),
    ("chandra-12", Groebner.Examples.chandran(12, k=k)),
    ("chandra-13", Groebner.Examples.chandran(13, k=k)),
    ("cyclic-7", Groebner.Examples.cyclicn(7, k=k)),
    ("cyclic-8", Groebner.Examples.cyclicn(8, k=k)),
    ("cyclic-9", Groebner.Examples.cyclicn(9, k=k)),
    ("eco-11", Groebner.Examples.econ(11, k=k)),
    ("eco-12", Groebner.Examples.econ(12, k=k)),
    ("eco-13", Groebner.Examples.econ(13, k=k)),
    ("eco-14", Groebner.Examples.econ(14, k=k)),
    ("noon-7", Groebner.Examples.noonn(7, k=k)),
    ("noon-8", Groebner.Examples.noonn(8, k=k)),
    ("noon-9", Groebner.Examples.noonn(9, k=k)),
    ("noon-10", Groebner.Examples.noonn(10, k=k)),
    ("henrion-6", Groebner.Examples.henrion6(k=k)),
    ("henrion-7", Groebner.Examples.henrion7(k=k)),
    ("henrion-8", Groebner.Examples.henrion8(k=k)),
    ("katsura-9", Groebner.Examples.katsuran(9, k=k)),
    ("katsura-10", Groebner.Examples.katsuran(10, k=k)),
    ("katsura-11", Groebner.Examples.katsuran(11, k=k)),
    ("katsura-12", Groebner.Examples.katsuran(12, k=k)),
    ("reimer-6", Groebner.Examples.reimern(6, k=k)),
    ("reimer-7", Groebner.Examples.reimern(7, k=k)),
    ("reimer-8", Groebner.Examples.reimern(8, k=k)),
    ("hexapod", Groebner.Examples.hexapod(k=k)),
]
sian = [
    ("cholera", Groebner.Examples.Cholera(k=k)),
    ("hiv2", Groebner.Examples.HIV2(k=k)),
    ("goodwin", Groebner.Examples.Goodwin_with_weights(k=k)),
    ("crn", Groebner.Examples.ChemicalReactionNetwork(k=k)),

]
other = [
    ("alea6", Groebner.Examples.alea6(k=k)),
    ("jason210", Groebner.Examples.jason210(k=k)),
    ("gametwo2", Groebner.Examples.gametwo2(k=k)),
    ("yang1", Groebner.Examples.yang1(k=k)),
    ("bayes148", Groebner.Examples.bayes148(k=k)),
    ("mayr42", Groebner.Examples.mayr42(k=k))
]
random = sort(reduce(vcat, [
    ("rand-$(n)-$d", randsys(n, d))
    for n in 4:6, d in 4:6
]), by=first)
random = sort(reduce(vcat, [
    ("rand-$(n)-$d", randsys(n, d))
    for n in 8:14, d in 2:2
]), by=first)
systems = vcat(system_solving, sian, other, random)
systems = vcat(random)

println("Running the following systems: ", map(first, systems))

timers = []
for (name, sys) in systems
    @info "Running $name.."
    # t1, t2, t3, t4, t5 = 0, 0, 0, 0, 0
    for t in 1:2
        TimerOutputs.enable_timer!(Groebner._TIMER); 
        reset_timer!(Groebner._TIMER); 
        groebner(sys); 
        if t == 2 show(Groebner._TIMER, allocations=false); println() end
    end
    push!(timers, [name, copy(Groebner._TIMER)])
end

labels = ["Name", "Select Pairs", "Symb Prepr", "Reduction", "Update", "Autoreduce", "IO", "Total"]
data = []
for entry in timers
    name = entry[1]
    timer = entry[2]
    time_select_pairs = TimerOutputs.time(timer["f4!"]["f4_select_critical_pairs!"])
    time_symb_preprc = TimerOutputs.time(timer["f4!"]["f4_symbolic_preprocessing!"])
    time_reduction = TimerOutputs.time(timer["f4!"]["f4_reduction!"])
    time_update = TimerOutputs.time(timer["f4!"]["f4_update!"])
    time_autoreduce = TimerOutputs.time(timer["f4!"]["f4_autoreduce!"])
    time_io = TimerOutputs.time(timer["io_convert_polynomials_to_ir"]) +
              TimerOutputs.time(timer["io_convert_ir_to_polynomials"]) +
              TimerOutputs.time(timer["ir_convert_internal_to_ir"]) +
                TimerOutputs.time(timer["ir_convert_ir_to_internal"])
    time_total = TimerOutputs.tottime(timer)
    push!(data, [name, time_select_pairs, time_symb_preprc, time_reduction, time_update, time_autoreduce, time_io, time_total])
end

matrix = permutedims(reduce(hcat, data))

threshold = 0.15
hl_select = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 2 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_symb_prepr = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 3 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_update = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 5 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_io = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 7 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"bold bg:light_gray")

pretty_table(
    matrix, 
    column_labels=labels,
    title="Timings in seconds",
    formatters=[
        (v,i,j) -> v isa Number ? @sprintf("%.2f", v / 1e9) : v, 
    ],
    row_group_labels = [1 => "System solving", length(system_solving)+1 => "SIAN", length(system_solving)+length(sian)+1 => "Other"],
    highlighters  = [hl_select, hl_symb_prepr, hl_update, hl_io],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    summary_row_labels = ["Total"],
    summary_rows = [(data, i) -> i > 1 ? @sprintf("%.2f", sum(data[:, i]) / 1e9) : ""],
    fit_table_in_display_vertically = false,
)
