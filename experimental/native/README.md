### Native codes!

1)

```julia
using Groebner, InteractiveUtils

io = open("llvm_vec4.txt", "w")
code_llvm(
    io,
    Groebner.reduce_dense_row_by_sparse_row!,
    Tuple{
        Vector{Groebner.CompositeInt{4, Int64}},
        Vector{Int32},
        Vector{Groebner.CompositeInt{4, Int32}},
        Groebner.SignedCompositeArithmeticZp{
            Groebner.CompositeInt{4, Int64},
            Groebner.CompositeInt{4, Int32},
            Int64,
            4
        },
    },
    debuginfo=:none
)       
close(io)
```

2)

```julia
using Groebner, InteractiveUtils

io = open("native_vec4.txt", "w")
code_native(
    io,
    Groebner.reduce_dense_row_by_sparse_row!,
    Tuple{
        Vector{Groebner.CompositeInt{4, Int64}},
        Vector{Int32},
        Vector{Groebner.CompositeInt{4, Int32}},
        Groebner.SignedCompositeArithmeticZp{
            Groebner.CompositeInt{4, Int64},
            Groebner.CompositeInt{4, Int32},
            Int64,
            4
        },
    },
    debuginfo=:none
)       
close(io)
```

3)

```julia
versioninfo()
```
