using BenchmarkTools

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

_setup1(n) = begin
    x = rand(UInt8.([0, 0, 0, 1, 2, 3]), n)
    y = rand(UInt8.([0, 0, 0, 1, 2, 3]), n)
    x, y
end
_setup2(n) = begin
    x = vcat(zeros(UInt8, n), rand(UInt8.([0, 0, 0, 1, 2, 3]), n))
    y = vcat(zeros(UInt8, n), rand(UInt8.([0, 0, 0, 1, 2, 3]), n))
    x, y
end
_setup3(T, n) = begin
    s = rand(T.([0, 0, 0, 1, 2, 3]), 3)
    x = Groebner.monom_construct_from_vector(
        Groebner.ExponentVector{T},
        vcat(zeros(T, n), s)
    )
    y = Groebner.monom_construct_from_vector(
        Groebner.ExponentVector{T},
        vcat(zeros(T, n), reverse(s))
    )
    z = similar(x)
    @assert Groebner.monom_totaldeg(x) == Groebner.monom_totaldeg(y)
    z, x, y
end

#########

function scalar_cmp_lex(x, y)
    i = 2
    @inbounds while i < length(x) && x[i] == y[i]
        i += 1
    end
    @inbounds x[i] < y[i]
end

begin
    n, step = 1, 5
    while n < 500
        @info "n = $n"
        print("scalar_cmp_lex\t\t\t")
        @btime scalar_cmp_lex(xx, yy) setup = begin
            cc, xx, yy = _setup3(Int8, max(1, $n))
        end
        print("Groebner.monom_is_equal\t\t")
        @btime Groebner.monom_is_equal(xx, yy) setup = begin
            cc, xx, yy = _setup3(Int8, max(1, $n))
        end
        print("Groebner.monom_copy\t\t")
        @btime Groebner.monom_copy(xx) setup = begin
            cc, xx, yy = _setup3(Int8, max(1, $n))
        end
        print("Groebner.monom_is_divisible\t")
        @btime Groebner.monom_is_divisible(xx, yy) setup = begin
            cc, xx, yy = _setup3(Int8, $n)
        end
        print("Groebner.monom_product!\t\t")
        @btime Groebner.monom_product!(cc, xx, yy) setup = begin
            cc, xx, yy = _setup3(Int8, max(1, $n))
        end
        print("Groebner.monom_lcm!\t\t")
        @btime Groebner.monom_lcm!(cc, xx, yy) setup = begin
            cc, xx, yy = _setup3(Int8, max(1, $n))
        end
        print("Groebner.monom_is_gcd_const\t")
        @btime Groebner.monom_is_gcd_const(xx, yy) setup = begin
            cc, xx, yy = _setup3(Int8, $n)
        end
        print("Groebner.monom_isless:lex\t")
        @btime Groebner.monom_isless(xx, yy, _ord) setup = begin
            cc, xx, yy = _setup3(Int8, max(1, $n))
            _ord = Groebner._Lex{true}(ones(Int, length(xx)))
        end
        print("Groebner.monom_isless:drl\t")
        @btime Groebner.monom_isless(xx, yy, _ord) setup = begin
            cc, xx, yy = map(reverse, _setup3(Int8, max(1, $n)))
            xx[1] = xx[end]
            yy[1] = yy[end]
            @assert Groebner.monom_totaldeg(xx) == Groebner.monom_totaldeg(yy)
            _ord = Groebner._DegRevLex{true}(ones(Int, length(xx)))
        end
        n += step
        step = ceil(Int, step * 1.2)
    end
end

_setup4(n) = begin
    s = rand(UInt8.([0, 0, 0, 1, 2, 3]), 2)
    a, b = vcat(zeros(UInt8, n), s), vcat(zeros(UInt8, n), reverse(s))
    a, b = reverse(a), reverse(b)
    x = Groebner.monom_construct_from_vector(Groebner.ExponentVector{UInt8}, a)
    y = Groebner.monom_construct_from_vector(Groebner.ExponentVector{UInt8}, b)
    vT(n) =
        if n + 2 < 8
            Groebner.PackedTuple1
        elseif n + 2 < 16
            Groebner.PackedTuple2
        elseif n + 2 < 24
            Groebner.PackedTuple3
        end
    xpacked = Groebner.monom_construct_from_vector(vT(n){UInt64, UInt8}, a)
    ypacked = Groebner.monom_construct_from_vector(vT(n){UInt64, UInt8}, b)
    x, y, xpacked, ypacked
end
begin
    n, step = 1, 3
    while n < 22
        @info "n = $n"
        print("Groebner.monom_isless:drl:packed\t")
        @btime Groebner.monom_isless(xpacked, ypacked, _ord) setup = begin
            x, y, xpacked, ypacked = _setup4($n)
            _ord = Groebner._DegRevLex{true}(ones(Int, length(x)))
            @assert Groebner.monom_totaldeg(xpacked) == Groebner.monom_totaldeg(ypacked)
            tmp1, tmp2 =
                Vector{Int8}(undef, length(x) - 1), Vector{Int8}(undef, length(x) - 1)
            @assert Groebner.monom_to_vector!(tmp1, x) ==
                    Groebner.monom_to_vector!(tmp2, xpacked)
            @assert Groebner.monom_to_vector!(tmp1, y) ==
                    Groebner.monom_to_vector!(tmp2, ypacked)
        end
        print("Groebner.monom_isless:drl:expvect\t")
        @btime Groebner.monom_isless(x, y, _ord) setup = begin
            x, y, xpacked, ypacked = _setup4($n)
            _ord = Groebner._DegRevLex{true}(ones(Int, length(x)))
            @assert Groebner.monom_totaldeg(x) == Groebner.monom_totaldeg(y)
        end
        n += step
    end
end

begin
    n, step = 1, 5
    while n < 500
        @info "n = $n"
        for _ in 1:1_000
            x, y = _setup3(n)
            res1 = vector_are_orth(x, y)
            res2 = _vec_check_orth(x, y)
            @assert res1 == res2
        end
        @btime vector_are_orth(xx, yy) setup = begin
            xx, yy = _setup3($n)
        end
        @btime _vec_check_orth(xx, yy) setup = begin
            xx, yy = _setup3($n)
        end
        n += step
        step = ceil(Int, step * 1.2)
    end
end

begin
    _setup1(n) = begin
        x = rand(UInt8.([0, 0, 0, 1, 2, 3]), n)
        y = rand(UInt8.([0, 0, 0, 1, 2, 3]), n)
        x, y
    end
    _setup2(n) = begin
        x = vcat(zeros(UInt8, n), rand(UInt8.([0, 0, 0, 1, 2, 3]), n))
        y = vcat(zeros(UInt8, n), rand(UInt8.([0, 0, 0, 1, 2, 3]), n))
        x, y
    end
    _setup3(n) = begin
        x = vcat(zeros(Int16, n), rand(Int16.([0, 0, 0, 1, 2, 3]), n))
        y = vcat(zeros(Int16, n), rand(Int16.([0, 0, 0, 1, 2, 3]), n))
        x, y
    end
    n, step = 1, 5
    while n < 500
        @info "n = $n"
        n += step
        step = ceil(Int, step * 1.2)
        for _ in 1:100
            x, y = _setup3(n)
            res1 = vector_any_lt(x, y)
            res2 = vector_any_lt_simd(x, y)
            @assert res1 == res2
        end
        @btime vector_any_lt(xx, yy) setup = begin
            xx, yy = _setup3($n)
        end
        @btime vector_any_lt_simd(xx, yy) setup = begin
            xx, yy = _setup3($n)
        end
    end
end

#########

begin
    _setup1(n) = begin
        x = rand(UInt8.([0, 0, 0, 1, 2, 3]), n)
        y = rand(UInt8.([0, 0, 0, 1, 2, 3]), n)
        similar(x), x, y
    end
    n, step = 1, 5
    while n < 500
        @info "n = $n"
        n += step
        step = ceil(Int, step * 1.2)
        for _ in 1:100
            c, x, y = _setup1(n)
            res1 = vector_emax_1!(copy(c), x, y)
            res2 = vector_emax_2!(copy(c), x, y)
            @assert res1 == res2
        end
        @btime vector_emax_1!(cc, xx, yy) setup = begin
            cc, xx, yy = _setup1($n)
        end
        @btime vector_emax_2!(cc, xx, yy) setup = begin
            cc, xx, yy = _setup1($n)
        end
    end
end
