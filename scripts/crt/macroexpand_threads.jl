#=
@macroexpand Base.Threads.@threads for i in 1:100
    println(i)
end
=#

# quote
#     local var"#17#threadsfor_fun"
#     begin
#         let var"#16#range" = 1:100
#             function var"#17#threadsfor_fun"(var"#28#tid" = 1; onethread = false)
#                 var"#21#r" = var"#16#range"
#                 var"#22#lenr" = Base.Threads.length(var"#21#r")
#                 if onethread
#                     var"#28#tid" = 1
#                     (var"#23#len", var"#24#rem") = (var"#22#lenr", 0)
#                 else
#                     (var"#23#len", var"#24#rem") = Base.Threads.divrem(var"#22#lenr", Base.Threads.threadpoolsize())
#                 end
#                 if var"#23#len" == 0
#                     if var"#28#tid" > var"#24#rem"
#                         return
#                     end
#                     (var"#23#len", var"#24#rem") = (1, 0)
#                 end
#                 var"#25#f" = Base.Threads.firstindex(var"#21#r") + (var"#28#tid" - 1) * var"#23#len"
#                 var"#26#l" = (var"#25#f" + var"#23#len") - 1
#                 if var"#24#rem" > 0
#                     #= threadingconstructs.jl:242 =#
#                     if var"#28#tid" <= var"#24#rem"
#                         #= threadingconstructs.jl:243 =#
#                         var"#25#f" = var"#25#f" + (var"#28#tid" - 1)
#                         #= threadingconstructs.jl:244 =#
#                         var"#26#l" = var"#26#l" + var"#28#tid"
#                     else
#                         #= threadingconstructs.jl:246 =#
#                         var"#25#f" = var"#25#f" + var"#24#rem"
#                         #= threadingconstructs.jl:247 =#
#                         var"#26#l" = var"#26#l" + var"#24#rem"
#                     end
#                 end
#                 #= threadingconstructs.jl:251 =#
#                 for var"#27#i" = var"#25#f":var"#26#l"
#                     #= threadingconstructs.jl:252 =#
#                     local i = begin
#                                 $(Expr(:inbounds, true))
#                                 local var"#29#val" = var"#21#r"[var"#27#i"]
#                                 $(Expr(:inbounds, :pop))
#                                 var"#29#val"
#                             end
#                     #= threadingconstructs.jl:253 =#
#                     begin
#                         #= REPL[3]:2 =#
#                         println(i)
#                         #= REPL[3]:3 =#
#                     end
#                     #= threadingconstructs.jl:254 =#
#                 end
#             end
#         end
#     end
#     Base.Threads.threading_run(var"#17#threadsfor_fun", false)
# end

using Base.Threads: nthreads, @threads, @spawn
using Base.Iterators: partition

tasks_per_thread = 2

some_data = collect(1:1000)

chunk_size = max(1, length(some_data) ÷ (tasks_per_thread * nthreads()))
data_chunks = partition(some_data, chunk_size) # partition your data into chunks that
                                               # individual tasks will deal with
#See also ChunkSplitters.jl and SplittablesBase.jl for partitioning data

tasks = Task[]
for chunk in data_chunks
    @info "chunk" chunk
    task = @spawn begin
        # state = some_initial_value
        for x in chunk
            # Core.print("threadid=$(Base.Threads.threadid()), value=$x\n")
            # state = some_operator(state, f(x))
        end
        # return state
    end
    push!(tasks, task)
end
fetch.(tasks)
