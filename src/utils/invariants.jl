# Custom assertions for Groebner
#
# Provides the @invariant macro

"""
    @invariant expr

Check that `expr` evaluates to `true` at runtime. If not, throw an
`AssertionError`.

It should be fine to use the `@invariant` macro in performance-critical code, as
long as Groebner is compiled with `invariants_enabled() = false`.

## Examples

```jldoctest
@invariant 2 > 1
@invariant 1 > 2  # throws!
```
"""
macro invariant(arg)
    # NOTE: This formulation of @invariant does not propagate the error source
    # path out of the macro
    expr = quote
        if $(@__MODULE__).invariants_enabled()
            @assert $arg
        else
            nothing
        end
    end
    esc(expr)
end
