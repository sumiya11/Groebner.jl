# Custom assertions for Groebner
#
# Provides the @invariant macro that can be turned off to have no runtime
# overhead.

"""
    @invariant expr

Check that `expr` evaluates to `true` at runtime. If not, throw an
`AssertionError`.

## Examples

```jldoctest
@invariant 2 > 1
@invariant 1 > 2  # throws!
```
"""
macro invariant(arg)
    esc(:(
        if $(@__MODULE__).invariants_enabled()
            @assert $arg
        else
            nothing
        end
    ))
    # NOTE: This formulation of @invariant does not propagate the error source
    # path out of the macro
    # if !invariants_enabled()
    #     return nothing
    # end
    # esc(:(@assert $arg))
end
