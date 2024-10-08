# This file is a part of Groebner.jl. License is GNU GPL v2.

# Provides the @invariant macro for custom assertion checking.

"""
    @invariant expr

Check that `expr` evaluates to `true` at runtime. If not, throw an
`AssertionError`.

It should be fine to use the `@invariant` macro in performance-critical code, as
long as the code is compiled with `Groebner.invariants_enabled() = false`.
"""
macro invariant(arg)
    expr = quote
        if $(@__MODULE__).invariants_enabled()
            @assert $arg
        else
            nothing
        end
    end
    esc(expr)
end

###
# DANGER # DANGER # DANGER #

"""
    unsafe_unreachable() --> Nothing

Invokes undefined behavior in LLVM.

The compiler may assume that this code is never reached during normal program
execution, and the program may (or may not) be optimized accordingly.

!!! danger
    Only use when able to prove that this is never called during program execution.

## Refs.

    - https://llvm.org/docs/LangRef.html#unreachable-instruction
    - https://llvm.org/docs/LangRef.html#llvm-assume-intrinsic
"""
unsafe_unreachable() = Base.llvmcall("unreachable", Cvoid, Tuple{})

"""
    unsafe_assume(condition) --> Nothing

May convince the compiler that the `condition` holds.

!!! danger
    Only use when able to prove that the condition always holds.
"""
function unsafe_assume(condition)
    # The idea is taken from UnsafeAssume.jl
    # https://gitlab.com/nsajko/UnsafeAssume.jl
    condition || unsafe_unreachable()
    nothing
end
