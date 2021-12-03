

# adds information about executed line to standard @time macro
# The code by Tamas Papp from
#   https://discourse.julialang.org/t/modifying-the-time-macro/2790
macro vtime(ex)
    quote
        local stats = Base.gc_num()
        local elapsedtime = time_ns()
        local val = $(esc(ex))
        elapsedtime = time_ns() - elapsedtime
        local diff = Base.GC_Diff(Base.gc_num(), stats)
        println($(Meta.quot(ex))) # print
        Base.time_print(elapsedtime, diff.allocd, diff.total_time,
                        Base.gc_alloc_count(diff))
        println()
        val
    end
end
