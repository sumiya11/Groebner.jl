### Native codes!

1)

```julia
using Groebner, InteractiveUtils

io = open("llvm_vec4.txt", "w")
code_llvm(
    io,
    Groebner.linalg_vector_addmul_sparsedense!,
    Tuple{
        Vector{Groebner.CompositeNumber{4, Int64}},
        Vector{Int32},
        Vector{Groebner.CompositeNumber{4, Int32}},
        Groebner.SignedCompositeArithmeticZp{
            Groebner.CompositeNumber{4, Int64},
            Groebner.CompositeNumber{4, Int32},
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
    Groebner.linalg_vector_addmul_sparsedense!,
    Tuple{
        Vector{Groebner.CompositeNumber{4, Int64}},
        Vector{Int32},
        Vector{Groebner.CompositeNumber{4, Int32}},
        Groebner.SignedCompositeArithmeticZp{
            Groebner.CompositeNumber{4, Int64},
            Groebner.CompositeNumber{4, Int32},
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
using InteractiveUtils; versioninfo()
```
