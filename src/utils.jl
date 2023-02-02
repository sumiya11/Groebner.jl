
# This function is adapted from HomotopyContinuation.jl 
"""
@tspawnat tid -> task
Mimics `Base.Threads.@spawn`, but assigns the task to thread `tid`.
# Example
```julia
julia> t = @tspawnat 4 Threads.threadid()
Task (runnable) @0x0000000010743c70
julia> fetch(t)
4
```
"""
macro tspawnat(thrdid, expr)
    thunk = esc(:(() -> ($expr)))
    var = esc(Base.sync_varname)
    tid = esc(thrdid)
    quote
        if $tid < 1 || $tid > Threads.nthreads()
            throw(
                AssertionError(
                    "@tspawnat thread assignment ($($tid)) must be between 1 and Threads.nthreads() (1:$(Threads.nthreads()))",
                ),
            )
        end
        local task = Task($thunk)
        task.sticky = false
        ccall(:jl_set_task_tid, Cvoid, (Any, Cint), task, $tid - 1)
        if $(Expr(:isdefined, var))
            push!($var, task)
        end
        schedule(task)
        task
    end
end
