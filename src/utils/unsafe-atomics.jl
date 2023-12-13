# Adapted from https://github.com/tkf/TSXPlayground.jl
# The license is MIT
module UnsafeAtomics

using Base.Sys: WORD_SIZE
using Base.Threads: inttypes, llvmtypes
using Core.Intrinsics: llvmcall

const unordered = Val{:unordered}()
const monotonic = Val{:monotonic}()
const acquire = Val{:acquire}()
const release = Val{:release}()
const acq_rel = Val{:acq_rel}()
const seq_cst = Val{:seq_cst}()

const orderings = [:unordered, :monotonic, :acquire, :release, :acq_rel, :seq_cst]

for typ in inttypes
    lt = llvmtypes[typ]
    rt = "$lt, $lt*"

    @eval @inline load(x::Ptr{$typ}) = load(x, seq_cst)
    for ord in orderings
        ord in [:release, :acq_rel] && continue

        @eval function load(x::Ptr{$typ}, ::$(Val{ord}))
            return llvmcall(
                $("""
                %ptr = inttoptr i$WORD_SIZE %0 to $lt*
                %rv = load atomic $rt %ptr $ord, align $(sizeof(typ))
                ret $lt %rv
                """),
                $typ,
                Tuple{Ptr{$typ}},
                x
            )
        end
    end

    @eval @inline store!(x::Ptr{$typ}, v::$typ) = store!(x, v, seq_cst)
    for ord in orderings
        ord in [:acquire, :acq_rel] && continue

        @eval function store!(x::Ptr{$typ}, v::$typ, ::$(Val{ord}))
            return llvmcall(
                $("""
                %ptr = inttoptr i$WORD_SIZE %0 to $lt*
                store atomic $lt %1, $lt* %ptr $ord, align $(sizeof(typ))
                ret void
                """),
                Cvoid,
                Tuple{Ptr{$typ}, $typ},
                x,
                v
            )
        end
    end

    @eval @inline cas!(x::Ptr{$typ}, cmp::$typ, new::$typ) =
        cas!(x, cmp, new, seq_cst, seq_cst)
    for success_ordering in orderings[2:end],
        failure_ordering in [:monotonic, :acquire, :seq_cst]

        @eval function cas!(
            x::Ptr{$typ},
            cmp::$typ,
            new::$typ,
            ::$(Val{success_ordering}),
            ::$(Val{failure_ordering})
        )
            return llvmcall(
                $(
                    """
                    %ptr = inttoptr i$WORD_SIZE %0 to $lt*
                    %rs = cmpxchg $lt* %ptr, $lt %1, $lt %2 $success_ordering $failure_ordering
                    %rv = extractvalue { $lt, i1 } %rs, 0
                    ret $lt %rv
                    """
                ),
                $typ,
                Tuple{Ptr{$typ}, $typ, $typ},
                x,
                cmp,
                new
            )
        end
    end

    for rmwop in [:add, :sub, :xchg, :and, :nand, :or, :xor, :max, :min]
        rmw = string(rmwop)
        fn = Symbol(rmw, "!")
        if (rmw == "max" || rmw == "min") && typ <: Unsigned
            # LLVM distinguishes signedness in the operation, not the integer type.
            rmw = "u" * rmw
        end
        @eval @inline $fn(x::Ptr{$typ}, v::$typ) = $fn(x, v, seq_cst)
        for ord in orderings
            @eval function $fn(x::Ptr{$typ}, v::$typ, ::$(Val{ord}))
                return llvmcall(
                    $("""
                    %ptr = inttoptr i$WORD_SIZE %0 to $lt*
                    %rv = atomicrmw $rmw $lt* %ptr, $lt %1 $ord
                    ret $lt %rv
                    """),
                    $typ,
                    Tuple{Ptr{$typ}, $typ},
                    x,
                    v
                )
            end
        end
    end
end

end  # module UnsafeAtomics
