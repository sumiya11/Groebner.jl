using Groebner, InteractiveUtils

code_native(
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
