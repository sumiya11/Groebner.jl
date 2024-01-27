using BenchmarkTools

function addmul!(
    row::Vector{UInt64},
    indices::Vector{Int64},
    coeffs::Vector{UInt32},
    mul::UInt64
)
    @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx] + mul * UInt64(coeffs[j])
        row[idx] = a
    end
    nothing
end

function addmul_2!(row, indices, coeffs, mul)
    textir = """
    ; row, indices, coeffs, mul, len
    define void @entry(i64 %0, i64 %1, i64 %2, i64 %3, i64 %4) #0 {
    top:
        %row = inttoptr i64 %0 to i64*
        %inds = inttoptr i64 %1 to i64*
        %coeffs = inttoptr i64 %2 to i32*
        %lenm7 = add nsw i64 %4, -3
        %dosimditer = icmp ugt i64 %4, 3
        br i1 %dosimditer, label %L9.lr.ph, label %L32

    L9.lr.ph:
        %len8 = and i64 %4, 9223372036854775804
        br label %L9

    L9:
        %i = phi i64 [ 0, %L9.lr.ph ], [ %vinc, %L30 ]
        %indsi = getelementptr inbounds i64, i64* %inds, i64 %i
        %idx = load i64, i64* %indsi, align 8
        %coeffsi = getelementptr inbounds i32, i32* %coeffs, i64 %i
        %coeff = load i32, i32* %coeffsi, align 8
        %idxm1 = add i64 %idx, -1
        %rowi = getelementptr inbounds i64, i64* %row, i64 %idxm1
        %x = load i64, i64* %rowi, align 8
        %coeff64 = zext i32 %coeff to i64
        %mulcoeff = mul i64 %coeff64, %3
        %newx = add i64 %x, %mulcoeff
        store i64 %newx, i64* %rowi, align 8
        ;
        %i2 = add i64 %i, 1
        %indsi2 = getelementptr inbounds i64, i64* %inds, i64 %i2
        %idx2 = load i64, i64* %indsi2, align 8
        %coeffsi2 = getelementptr inbounds i32, i32* %coeffs, i64 %i2
        %coeff2 = load i32, i32* %coeffsi2, align 8
        %idxm12 = add i64 %idx2, -1
        %rowi2 = getelementptr inbounds i64, i64* %row, i64 %idxm12
        %x2 = load i64, i64* %rowi2, align 8
        %coeff642 = zext i32 %coeff2 to i64
        %mulcoeff2 = mul i64 %coeff642, %3
        %newx2 = add i64 %x2, %mulcoeff2
        store i64 %newx2, i64* %rowi2, align 8
        ;
        %i3 = add i64 %i2, 1
        %indsi3 = getelementptr inbounds i64, i64* %inds, i64 %i3
        %idx3 = load i64, i64* %indsi3, align 8
        %coeffsi3 = getelementptr inbounds i32, i32* %coeffs, i64 %i3
        %coeff3 = load i32, i32* %coeffsi3, align 8
        %idxm13 = add i64 %idx3, -1
        %rowi3 = getelementptr inbounds i64, i64* %row, i64 %idxm13
        %x3 = load i64, i64* %rowi3, align 8
        %coeff643 = zext i32 %coeff3 to i64
        %mulcoeff3 = mul i64 %coeff643, %3
        %newx3 = add i64 %x3, %mulcoeff3
        store i64 %newx3, i64* %rowi3, align 8
        ;
        %i4 = add i64 %i3, 1
        %indsi4 = getelementptr inbounds i64, i64* %inds, i64 %i4
        %idx4 = load i64, i64* %indsi4, align 8
        %coeffsi4 = getelementptr inbounds i32, i32* %coeffs, i64 %i4
        %coeff4 = load i32, i32* %coeffsi4, align 8
        %idxm14 = add i64 %idx4, -1
        %rowi4 = getelementptr inbounds i64, i64* %row, i64 %idxm14
        %x4 = load i64, i64* %rowi4, align 8
        %coeff644 = zext i32 %coeff4 to i64
        %mulcoeff4 = mul i64 %coeff644, %3
        %newx4 = add i64 %x4, %mulcoeff4
        store i64 %newx4, i64* %rowi4, align 8
        ;
        br label %L30

    L30:
        %vinc = add nuw nsw i64 %i, 4
        %continue = icmp slt i64 %vinc, %lenm7
        br i1 %continue, label %L9, label %L32

    L32:
        %cumi = phi i64 [ 0, %top ], [ %len8, %L30 ]
        %done = icmp eq i64 %cumi, %4
        br i1 %done, label %common.ret, label %L51

    common.ret:
        ret void

    L51:
        %si = phi i64 [ %inc, %L67 ], [ %cumi, %L32 ]
        %sindsi = getelementptr inbounds i64, i64* %inds, i64 %si
        %sidx = load i64, i64* %sindsi, align 8
        %scoeffsi = getelementptr inbounds i32, i32* %coeffs, i64 %si
        %scoeff = load i32, i32* %scoeffsi, align 8
        %sidxm1 = add i64 %sidx, -1
        %srowi = getelementptr inbounds i64, i64* %row, i64 %sidxm1
        %sx = load i64, i64* %srowi, align 8
        %scoeff64 = zext i32 %scoeff to i64
        %smulcoeff = mul i64 %scoeff64, %3
        %snewx = add i64 %sx, %smulcoeff
        store i64 %snewx, i64* %srowi, align 8
        br label %L67

    L67:
        %inc = add i64 %si, 1
        %dobreak = icmp eq i64 %inc, %4
        br i1 %dobreak, label %common.ret, label %L51

    }
    attributes #0 = { alwaysinline }
"""

    GC.@preserve row indices coeffs begin
        Base.llvmcall(
            (textir, "entry"),
            Nothing,
            Tuple{Ptr{UInt64}, Ptr{Int64}, Ptr{UInt32}, UInt64, Int64},
            Base.unsafe_convert(Ptr{UInt64}, row),
            Base.unsafe_convert(Ptr{Int64}, indices),
            Base.unsafe_convert(Ptr{UInt32}, coeffs),
            mul,
            length(indices)
        )
    end
end

@code_native debuginfo = :none addmul!(UInt[], Int[], UInt32[], UInt(1))
@code_native debuginfo = :none addmul_2!(UInt[], Int[], UInt32[], UInt(1))

n = 10_000
for k in 1:5:200
    @info "n = $n, k = $k"
    for _ in 1:100
        row = rand(UInt.(0:100), n)
        indices = sort(unique(rand(1:(n), k)))
        coeffs = rand(UInt32.(1:100), length(indices))
        mul = rand(UInt)
        crow = copy(row)
        # println(row[indices])
        addmul!(row, indices, coeffs, mul)
        addmul_2!(crow, indices, coeffs, mul)
        # println(row[indices])
        # println(crow[indices])
        # println(indices)
        # println(coeffs)
        @assert row == crow
    end
    @btime addmul!(row, indices, coeffs, mul) setup = begin
        row = rand(UInt, $n)
        indices = sort(unique(rand(1:($n), $k)))
        coeffs = rand(UInt32, length(indices))
        mul = rand(UInt)
    end
    @btime addmul_2!(row, indices, coeffs, mul) setup = begin
        row = rand(UInt, $n)
        indices = sort(unique(rand(1:($n), $k)))
        coeffs = rand(UInt32, length(indices))
        mul = rand(UInt)
    end
end
