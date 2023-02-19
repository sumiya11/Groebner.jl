include((@__DIR__) * "../../benchmark/systems/MQ/parser.jl")

system = read_MQ_GF("mq_n15_m30_p2_s0")

VSCodeServer.@profview gb = Groebner.groebner(system, linalg=:prob);
