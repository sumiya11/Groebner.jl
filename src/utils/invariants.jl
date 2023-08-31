# Custom assertions for Groebner
#
# Provides the @invariant macro

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
    # NOTE: This formulation of @invariant does not propagate the error source
    # path out of the macro
    esc(:(
        if $(@__MODULE__).invariants_enabled()
            @assert $arg
        else
            nothing
        end
    ))
    # if !invariants_enabled()
    #     return nothing
    # end
    # esc(:(@assert $arg))
end
