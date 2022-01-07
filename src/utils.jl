

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
        Logging.# @debug $(Meta.quot(ex)) # print
        Base.time_print(elapsedtime, diff.allocd, diff.total_time,
                        Base.gc_alloc_count(diff))
        println()
        val
    end
end


# wraps nested loops into one for convenience
#   @uwu for (i, j) in (arr1,arr2)
# is equivalent to
#   for i in arr1
#       for j in arr2
macro uwu(forloop)
    head, body = forloop.args
    iters, iterables = head.args

    stri   = map(string, iters.args)

    strhead = reduce(*,
        map(x -> "for "*x[1]*" in "*string(:($(x[2])))*" ",
        zip(stri, iterables.args))
    )

    strbody = string(body)
    whole = strhead*"\n"*strbody*"\nend"^length(stri)

    Meta.parse(whole)
end
