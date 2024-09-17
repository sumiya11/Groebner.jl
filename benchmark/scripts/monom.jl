using BenchmarkTools, PrettyTables, IOCapture

using HostCPUFeatures
using HostCPUFeatures:
    register_size,
    pick_vector_width,
    pick_vector_width_shift,
    simd_integer_register_size,
    fma_fast,
    has_feature,
    register_count,
    cpu_name,
    register_size

#########

function setup_data(n, T, kind=:randeq)
    if kind === :randeq
        a, b = zeros(T, n), zeros(T, n)
        inds = rand(1:n, 10)
        a[inds] .= 1
        b[inds] .= 1
    end
    a, b
end

function setup_dense(n, T)
    a, b = setup_data(n, T)
    x = Groebner.monom_construct_from_vector(Groebner.ExponentVector{T}, a)
    y = Groebner.monom_construct_from_vector(Groebner.ExponentVector{T}, b)
    x, y
end

packed_type(n) =
    if n < 8
        Groebner.PackedTuple1
    elseif n < 16
        Groebner.PackedTuple2
    elseif n < 24
        Groebner.PackedTuple3
    elseif n < 32
        Groebner.PackedTuple4
    else
        nothing
    end

function setup_packed(n, T)
    a, b = setup_data(n, T)
    type = packed_type(n)
    isnothing(type) && return nothing
    x_packed = Groebner.monom_construct_from_vector(type{UInt64, T}, a)
    y_packed = Groebner.monom_construct_from_vector(type{UInt64, T}, b)
    x_packed, y_packed
end

function setup_sparse(n, T)
    a, b = setup_data(n, T)
    x_sparse = Groebner.monom_construct_from_vector(Groebner.SparseExponentVector{T}, a)
    y_sparse = Groebner.monom_construct_from_vector(Groebner.SparseExponentVector{T}, b)
    x_sparse, y_sparse
end

funcs = ["monom_lcm!", "monom_is_divisible"]
execs = [
    ((n, setup) -> @btime Groebner.monom_lcm!(z, x, y) setup = begin
        x, y = $(setup)($n)
        z = Groebner.monom_copy(x)
    end),
    ((n, setup) -> @btime Groebner.monom_is_divisible(x, y) setup = begin
        x, y = $(setup)($n)
    end)
]
impls = ["dense:u8", "dense:u32", "packed", "sparse"]
setup = [
    n -> setup_dense(n, UInt8),
    n -> setup_dense(n, UInt32),
    n -> setup_packed(n, UInt8),
    n -> setup_sparse(n, UInt8)
]
ns = [7, 8, 15, 16, 23, 24, 31, 32, 63, 64, 127, 128, 255, 256]
# ns = [255, 256, 400, 511, 512, 1023, 1024]

begin
    for (k, func) in enumerate(funcs[2:2])
        table = Matrix{Any}(undef, (length(ns), length(impls)))
        exec = execs[k]
        @info "func = $func"
        for (i, n) in enumerate(ns)
            @info "n = $n"
            row = Vector{Any}(undef, length(impls))
            for (j, impl) in enumerate(impls)
                print("$func:$impl\t")
                if (setup[j])(n) === nothing
                    println("-")
                    row[j] = "-"
                    continue
                end
                c = IOCapture.capture() do
                    exec(n, setup[j])
                end
                s = strip(c.output)
                println(s)
                row[j] = parse(Float64, s[1:findfirst(' ', s)])
            end
            table[i, :] .= row
        end
        pretty_table(table, title=func, header=impls, tf=tf_markdown, row_labels=ns)
    end
end
