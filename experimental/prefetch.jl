
# Adapted from
#   https://discourse.julialang.org/t/how-to-measure-cache-misses-in-julia/20875/2

"""
    perf(f, args)

Run `perf` program with `args` while executing `f`.

# Examples
```
julia> A = randn(100, 100); B = randn(100, 100);

julia> perf(`stat`) do
           @benchmark \$A * \$B
       end
```
"""
function perf(f, args)
    pid = getpid()
    cmd = `perf $args --pid=$pid`
    proc = run(pipeline(cmd, stdout=stdout, stderr=stderr); wait=false)
    try
        return f()
    finally
        flush(stdout)
        flush(stderr)
        kill(proc, Base.SIGINT)
        wait(proc)
    end
end

perf_stat_l1(f, args = ``) =
    perf(f, `stat L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,L1-icache-load-misses $args`)

"""
    unsafe_prefetch(
        address::Union{Ptr,Integer},
        ::Val{rw} = Val(:read),
        ::Val{locality} = Val(0),
        ::Val{cache_type} = Val(:data),
    )
# Arguments
- `rw`: `:read` (0) or `:write` (1)
- `locality`: no locality/NTA (0) -- extremely local/T0 (3)
- `cache_type`: `:data` (1) or `:instruction` (0)
"""
unsafe_prefetch

@generated function unsafe_prefetch(
    address::Union{Ptr,Integer},
    ::Val{rw} = Val(:read),
    ::Val{locality} = Val(0),
    ::Val{cache_type} = Val(:data),
) where {locality,rw,cache_type}

    rw = get(Dict(:read => 0, :write => 1), rw, rw)
    cache_type = get(Dict(:data => 1, :instruction => 0), cache_type, cache_type)

    @assert rw in (0, 1)
    @assert locality in 0:3
    @assert cache_type in (0, 1)

    declaration = "declare void @llvm.prefetch(i8*, i32, i32, i32)"

    typ = (Int === Int64 ? "i64" : "i32")
    instructions = """
    %addr = inttoptr $typ %0 to i8*
    call void @llvm.prefetch(i8* %addr, i32 $rw, i32 $locality, i32 $cache_type)
    ret void
    """
    if VERSION < v"1.6.0-DEV.674"
        IR = (declaration, instructions)
    else
        IR = (
            """
            $declaration
            define void @entry($typ) #0 {
            top:
                $instructions
            }
            attributes #0 = { alwaysinline }
            """,
            "entry",
        )
    end

    quote
        $(Expr(:meta, :inline))
        Base.llvmcall($IR, Cvoid, Tuple{Ptr{Cvoid}}, address)
    end
end

@inline function tryprefetch(x::T) where {T}
    R = Core.Compiler.return_type(pointer, Tuple{T})
    if R !== Union{} && R <: Ptr
        unsafe_prefetch(pointer(x))
    end
end

