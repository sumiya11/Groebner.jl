# This file is a part of Groebner.jl. License is GNU GPL v2.

"""
    groebner(polynomials; options...)

Computes a Groebner basis of the ideal generated by `polynomials`.

## Arguments

- `polynomials`: an array of polynomials. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl.

## Returns

- `basis`: an array of polynomials, a Groebner basis.

## Possible Options

- `reduced`: A bool, if the returned basis must be autoreduced and unique.
  Default is `true`.
- `ordering`: Specifies the monomial ordering. Available monomial orderings are: 
    - `InputOrdering()` for inferring the ordering from the given `polynomials`
      (default), 
    - `Lex(args...)` for lexicographic, 
    - `DegLex(args...)` for degree lexicographic, 
    - `DegRevLex(args...)` for degree reverse lexicographic, 
    - `WeightedOrdering(args...)` for weighted ordering, 
    - `ProductOrdering(args...)` for block ordering, 
    - `MatrixOrdering(args...)` for matrix ordering.
    For details and examples see the corresponding documentation page.
- `certify`: A bool, whether to certify the obtained basis. When this option is
    `false`, the algorithm is randomized and the result is correct with high
    probability. When this option is `true`, the result is guaranteed to be
    correct in case the ideal is homogeneous. Default is `false`. 
- `linalg`: A symbol, linear algebra backend. Available options are: 
    - `:auto` for the automatic choice (default),
    - `:deterministic` for deterministic sparse linear algebra, 
    - `:randomized` for probabilistic sparse linear algebra.
- `threaded`: The use of multi-threading. Available options are: 
    - `:auto` for the automatic choice (default),
    - `:no` never use multi-threading, 
    - `:yes` allow the use of multi-threading.
    Additionally, it is possible to set the environment variable
    `GROEBNER_NO_THREADED` to `1` to disable all multi-threading in Groebner.jl.
    In this case, the environment variable takes precedence over the `threaded`
    option.
- `monoms`: Monomial representation used in the computations. Available options
  are: 
    - `:auto` for the automatic choice (default), 
    - `:dense` for classic dense exponent vectors,
    - `:packed` for packed representation. 
- `modular`: Modular computation algorithm. Only has effect when computing basis
    over rational numbers. Available options are:
    - `:auto` for the automatic choice (default),
    - `:classic_modular` for the classic multi-modular algorithm,
    - `:learn_and_apply` for the learn & apply algorithm.
- `seed`: The seed for randomization. Default is `42`.
- `homogenize`: Controls the use of homogenization in the algorithm. Available
  options are:
    - `:auto`, for the automatic choice (default).
    - `:yes`, always homogenize the input ideal,
    - `:no`, never homogenize the input ideal,

## Example

Using DynamicPolynomials.jl:

```jldoctest_
using Groebner, DynamicPolynomials
@polyvar x y
groebner([x*y^2 + x, y*x^2 + y])
```

Using AbstractAlgebra.jl:

```jldoctest_
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"]
groebner([x*y^2 + x, y*x^2 + y])
```

Using Nemo.jl:

```jldoctest_
using Groebner, Nemo
R, (x, y) = GF(2^30+3)["x", "y"]
groebner([x*y^2 + x, y*x^2 + y])
```

Or, say, in another monomial ordering:

```jldoctest_
# lex with y > x
groebner([x*y^2 + x, y*x^2 + y], ordering=Lex(y, x))

# degree reverse lexicographic
groebner([x*y^2 + x, y*x^2 + y], ordering=DegRevLex())
```

## Notes

- The function is thread-safe.
- For AbstractAlgebra.jl and Nemo.jl, the function is most efficient for
    polynomials over `GF(p)`, `Native.GF(p)`, and `QQ`.
- The default algorithm is probabilistic (with `certify=false`). Results are
    correct with high probability, however, no precise bound on the probability
    is known.
"""
function groebner(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:groebner, options)
    result = groebner0(polynomials, keywords)
    result
end

#=
    groebner(ring, monoms, coeffs; options...)

Computes a Groebner basis. This function is a part of low level interface.

## Arguments

- `ring`: polynomial ring, a `PolyRing` object.
- `monoms`: vector of monomials of input polynomials.
- `coeffs`: vector of coefficients of input polynomials.

## Returns

Returns a tuple (`gb_monoms`, `gb_coeffs`).

- `gb_monoms`: an array, the monomials of a Groebner basis.
- `gb_coeffs`: an array, the coefficients of a Groebner basis.

## Possible Options

Same as for public `groebner`.

## Example

Computing a Groebner basis of a bivariate ideal modulo 65537 in degrevlex:

```jldoctest_
using Groebner
ring = Groebner.PolyRing(2, Groebner.DegRevLex(), 65537)
# {x y - 1,  x^3 + 7 y^2}
monoms = [ [[1, 1], [0, 0]], [[3, 0], [0, 2]] ]
coeffs = [ [1,65536], [1, 7] ]
Groebner.groebner(ring, monoms, coeffs)
```

## Notes

Same as for `groebner`.
=#
function groebner(ring::PolyRing, monoms::AbstractVector, coeffs::AbstractVector; options...)
    keywords = KeywordArguments(:groebner, options)
    ir_is_valid_basic(ring, monoms, coeffs)
    result = groebner1(ring, monoms, coeffs, keywords)
    result
end

"""
    groebner_with_change_matrix(polynomials; options...)

Computes a Groebner basis of the ideal generated by `polynomials` and emits a
change matrix, that is, a map from the original generators to basis elements.

## Arguments

- `polynomials`: an array of polynomials. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl. For
    AbstractAlgebra.jl and Nemo.jl, coefficients of polynomials must belong to
    `GF(p)`, `Native.GF(p)`, or `QQ`. 

## Returns

Returns a tuple (`basis`, `matrix`).

- `basis`: an array of polynomials, a Groebner basis.
- `matrix`: a matrix, so that `matrix * polynomials == basis`.

## Possible Options

Same as for `groebner`.

## Example

Using AbstractAlgebra.jl:

```jldoctest_
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"]
f = [x*y^2 + x, y*x^2 + y]

g, m = groebner_with_change_matrix(f, ordering=DegRevLex())

@assert isgroebner(g, ordering=DegRevLex())
@assert m * f == g
```

## Notes

- Only `DegRevLex` ordering is supported.
- The function is thread-safe.
- The default algorithm is probabilistic (with `certify=false`). Results are
    correct with high probability, however, no precise bound on the probability
    is known.
"""
function groebner_with_change_matrix(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:groebner_with_change_matrix, options)
    result = groebner_with_change_matrix0(polynomials, keywords)
    result
end

"""
    groebner_learn(polynomials; options...)

Computes a Groebner basis of the ideal generated by `polynomials` and emits a
trace.

The trace can be used to speed up the computation of Groebner bases of
specializations of the same ideal as the one `groebner_learn` had been applied
to.

See also `groebner_apply!`.

## Arguments

- `polynomials`: an array of polynomials. Must be polynomials from
  AbstractAlgebra.jl or Nemo.jl over `GF(p)` or `Native.GF(p)`.

## Returns

Returns a tuple (`trace`, `basis`).

- `trace`: an object, a trace. Can be used in `groebner_apply!`.
- `basis`: an array of polynomials, a Groebner basis.

## Possible Options

Same as for `groebner`.

## Example

Using `groebner_learn` and `groebner_apply!` over the same ground field:

```jldoctest_
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Apply (same support, different coefficients)
flag, gb_2 = groebner_apply!(trace, [2x*y^2 + 3x, 4y*x^2 + 5y])

@assert flag
```

Using `groebner_learn` and `groebner_apply!` over different ground fields:

```jldoctest_
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y], ordering=DegRevLex())

# Create a ring with a different modulo
R2, (x2, y2) = GF(2^30+3)["x", "y"]

# Apply (different modulo)
flag, gb_2 = groebner_apply!(
    trace, 
    [2x2*y2^2 + 3x2, 4y2*x2^2 + 5y2], 
    ordering=DegRevLex()
)

@assert flag
@assert gb_2 == groebner([2x2*y2^2 + 3x2, 4y2*x2^2 + 5y2], ordering=DegRevLex())
```

Using `groebner_apply!` in batches:

```jldoctest_
using Groebner, AbstractAlgebra
R, (x, y) = polynomial_ring(GF(2^31-1), ["x", "y"], internal_ordering=:degrevlex)

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Create rings with some other moduli
R2, (x2, y2) = polynomial_ring(GF(2^30+3), ["x", "y"], internal_ordering=:degrevlex)
R3, (x3, y3) = polynomial_ring(GF(2^27+29), ["x", "y"], internal_ordering=:degrevlex)

# Two specializations of the same ideal
batch = ([2x2*y2^2 + 3x2, 4y2*x2^2 + 5y2], [4x3*y3^2 + 4x3, 5y3*x3^2 + 7y3])

# Apply for two sets of polynomials at once
flag, (gb_2, gb_3) = groebner_apply!(trace, batch)

@assert flag
@assert (gb_2, gb_3) == map(groebner, batch)
```

Perhaps, in a more involved example, we will compute Groebner bases of the
Katsura-9 system:

```jldoctest_
using Groebner, AbstractAlgebra, BenchmarkTools

# Create the system
kat = Groebner.Examples.katsuran(9, k=ZZ, internal_ordering=:degrevlex)

# Reduce the coefficients modulo 5 different primes
kat_0 = map(f -> map_coefficients(c -> GF(2^30 + 3)(c), f), kat)
kat_1 = map(f -> map_coefficients(c -> GF(2^30 + 7)(c), f), kat)
kat_2 = map(f -> map_coefficients(c -> GF(2^30 + 9)(c), f), kat)
kat_3 = map(f -> map_coefficients(c -> GF(2^30 + 15)(c), f), kat)
kat_4 = map(f -> map_coefficients(c -> GF(2^30 + 19)(c), f), kat)

# Learn the trace
trace, gb_0 = groebner_learn(kat_0);

# Compare the performance of applying with 1 input and with 4 different inputs:

# Apply for one system
@btime groebner_apply!(\$trace, \$kat_1);
#  46.824 ms (19260 allocations: 24.48 MiB)

# Apply for a batch of four systems
@btime groebner_apply!(\$trace, \$(kat_1, kat_2, kat_3, kat_4));
#  72.813 ms (23722 allocations: 59.44 MiB)
```

Observe the better amortized performance of the composite `groebner_apply!`.

## Notes

- The function is thread-safe.
"""
function groebner_learn(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:groebner_learn, options)
    result = groebner_learn0(polynomials, keywords)
    result
end

# This function is a part of low level interface.
function groebner_learn(ring::PolyRing, monoms::AbstractVector, coeffs::AbstractVector; options...)
    keywords = KeywordArguments(:groebner_learn, options)
    ir_is_valid_basic(ring, monoms, coeffs)
    result = groebner_learn1(ring, monoms, coeffs, keywords)
    result
end

"""
    groebner_apply!(trace, polynomials; options...)
    groebner_apply!(trace, batch::NTuple{N, Vector}; options...)

Computes a Groebner basis of the ideal generated by `polynomials` following the
given `trace`. 

See also `groebner_learn`.

## Arguments

- `trace`: a trace produced by `groebner_learn`.
- `polynomials`: an array of polynomials. Must be polynomials from
    AbstractAlgebra.jl or Nemo.jl over `GF(p)` or `Nemo.GF(p)`. It is possible
    to supply a tuple of `N` arrays of polynomials to compute `N` Groebner bases
    simultaneously. This could be more efficient overall than computing them in
    separate.

## Returns

Returns a tuple (`success`, `basis`).

- `success`: a bool, whether the call succeeded.
- `basis`: an array of polynomials, a Groebner basis.

## Possible Options

The `groebner_apply!` function automatically inherits most parameters from the
given `trace`.

## Example

For examples, see the documentation of `groebner_learn`.

## Notes

- In general, `success` may be a false positive. The probability of a false
  positive is considered to be low enough in some practical applications.

- This function is **not** thread-safe since it mutates `trace`.

- This function is **not** safe against program interruptions. For example,
    pressing `Ctrl + C` while `groebner_apply!(trace, ...)` is running may leave
    `trace` corrupted.
"""
function groebner_apply! end

# Specialization for a single input
function groebner_apply!(trace::WrappedTrace, polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:groebner_apply!, options)
    result = groebner_apply0!(trace, polynomials, keywords)
    result
end

const _supported_batch_size = (1, 2, 4, 8, 16, 32, 64, 128)

# Specialization for a batch of inputs
function groebner_apply!(
    trace::WrappedTrace,
    batch::NTuple{N, T}; # deliberately not ::Tuple{T, Vararg{T, Nminus1}}
    options...
) where {N, T <: AbstractVector}
    !(N in _supported_batch_size) &&
        throw(DomainError("The batch size must be one of the following: $_supported_batch_size"))
    keywords = KeywordArguments(:groebner_apply!, options)
    result = groebner_apply_batch0!(trace, batch, keywords)
    result
end

# Low level specialization for a single input
function groebner_apply!(
    trace::WrappedTrace,
    ring::PolyRing,
    monoms::AbstractVector,
    coeffs::AbstractVector;
    options...
)
    keywords = KeywordArguments(:groebner_apply!, options)
    ir_is_valid_basic(ring, monoms, coeffs)
    result = groebner_apply1!(trace, ring, monoms, coeffs, keywords)
    result
end

# Low level specialization for a batch of inputs
function groebner_apply!(trace::WrappedTrace, batch::NTuple{N, T}; options...) where {N, T}
    !(N in _supported_batch_size) &&
        throw(DomainError("The batch size must be one of the following: $_supported_batch_size"))
    keywords = KeywordArguments(:groebner_apply!, options)
    ir_is_valid_basic(batch)
    result = groebner_apply_batch1!(trace, batch, keywords)
    result
end

"""
    isgroebner(polynomials; options...)

Checks if `polynomials` forms a Groebner basis.

## Arguments

- `polynomials`: an array of polynomials. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl. 

## Returns

- `flag`: a bool, whether `polynomials` is a Groebner basis of the ideal
  generated by `polynomials`.

## Possible Options

- `ordering`: Specifies the monomial ordering. Available monomial orderings are: 
    - `InputOrdering()` for inferring the ordering from the given `polynomials`
      (default), 
    - `Lex()` for lexicographic, 
    - `DegLex()` for degree lexicographic, 
    - `DegRevLex()` for degree reverse lexicographic, 
    - `WeightedOrdering(weights)` for weighted ordering, 
    - `ProductOrdering(args...)` for block ordering, 
    - `MatrixOrdering(matrix)` for matrix ordering. 
  For details and examples see the corresponding documentation page.
- `certify`: a bool, whether to use a deterministic algorithm. Default is
  `false`.
- `seed`: The seed for randomization. Default value is `42`.

## Example

Using `DynamicPolynomials`:

```jldoctest_
using Groebner, DynamicPolynomials
@polyvar x y;
isgroebner([x*y^2 + x, y*x^2 + y])
```

Using `AbstractAlgebra`:

```jldoctest_
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"]
isgroebner([x*y^2 + x, y*x^2 + y])
```

## Notes

- The function is thread-safe.
- For AbstractAlgebra.jl and Nemo.jl, the function is most efficient for
    polynomials over `GF(p)`, `Native.GF(p)`, and `QQ`.
- The default algorithm is probabilistic (with `certify=false`). Results are
    correct with high probability, however, no precise bound on the probability
    is known. 
"""
function isgroebner(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:isgroebner, options)
    result = isgroebner0(polynomials, keywords)::Bool
    result
end

# This function is a part of low level interface.
function isgroebner(ring::PolyRing, monoms::AbstractVector, coeffs::AbstractVector; options...)
    keywords = KeywordArguments(:isgroebner, options)
    ir_is_valid_basic(ring, monoms, coeffs)
    result = isgroebner1(ring, monoms, coeffs, keywords)
    result
end

"""
    normalform(basis, to_be_reduced; options...)

Computes the normal form of polynomials `to_be_reduced` with respect to a
Groebner basis `basis`.

## Arguments

- `basis`: an array of polynomials, a Groebner basis. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl.
- `to_be_reduced`: either a single polynomial or an array of polynomials.
    Supports polynomials from AbstractAlgebra.jl, Nemo.jl, and
    DynamicPolynomials.jl.

## Returns

- `reduced`: either a single polynomial or an array of polynomials, the normal
  forms.

## Possible Options

- `check`: Check if `basis` forms a Groebner basis. Default is `true`.
- `ordering`: Specifies the monomial ordering. Available monomial orderings are: 
    - `InputOrdering()` for inferring the ordering from the given `polynomials`
      (default), 
    - `Lex()` for lexicographic, 
    - `DegLex()` for degree lexicographic, 
    - `DegRevLex()` for degree reverse lexicographic, 
    - `WeightedOrdering(weights)` for weighted ordering, 
    - `ProductOrdering(args...)` for block ordering, 
    - `MatrixOrdering(matrix)` for matrix ordering. 
    For details and examples see the corresponding documentation page.

## Example

Fining the normal form a single polynomial:

```jldoctest_
using Groebner, DynamicPolynomials
@polyvar x y;
normalform([y^2 + x, x^2 + y], x^2 + y^2 + 1)
```

Or, reducing two polynomials at a time:

```jldoctest_
using Groebner, DynamicPolynomials
@polyvar x y;
normalform([y^2 + x, x^2 + y], [x^2 + y^2 + 1, x^10*y^10])
```

## Notes

- The function is thread-safe.
- For AbstractAlgebra.jl and Nemo.jl, the function is most efficient for
    polynomials over `GF(p)`, `Native.GF(p)`, and `QQ`.
- The default algorithm is probabilistic (with `certify=false`). Results are
    correct with high probability, however, no precise bound on the probability
    is known.
"""
function normalform(basis::AbstractVector, to_be_reduced::AbstractVector; options...)
    keywords = KeywordArguments(:normalform, options)
    result = normalform0(basis, to_be_reduced, keywords)
    result
end

normalform(basis::AbstractVector, to_be_reduced; options...) =
    first(normalform(basis, [to_be_reduced]; options...))

function normalform(
    ring::PolyRing,
    monoms::AbstractVector,
    coeffs::AbstractVector,
    ring_to_be_reduced::PolyRing,
    monoms_to_be_reduced::AbstractVector,
    coeffs_to_be_reduced::AbstractVector;
    options...
)
    keywords = KeywordArguments(:normalform, options)
    ir_is_valid_basic(ring, monoms, coeffs)
    ir_is_valid_basic(ring_to_be_reduced, monoms_to_be_reduced, coeffs_to_be_reduced)
    result = normalform1(
        ring,
        monoms,
        coeffs,
        ring_to_be_reduced,
        monoms_to_be_reduced,
        coeffs_to_be_reduced,
        keywords
    )
    result
end

"""
    leading_term(polynomial; options...)

Returns the leading term of a polynomial.

## Arguments

- `polynomial`: a polynomial. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl.

## Possible Options

- `ordering`: Specifies the monomial ordering. Available monomial orderings are:
    - `InputOrdering()` for inferring the ordering from the given `polynomials`
      (default),
    - `Lex()` for lexicographic,
    - `DegLex()` for degree lexicographic,
    - `DegRevLex()` for degree reverse lexicographic,
    - `WeightedOrdering(weights)` for weighted ordering,
    - `ProductOrdering(args...)` for block ordering,
    - `MatrixOrdering(matrix)` for matrix ordering.
    For details and examples see the corresponding documentation page.

## Notes

- The function is thread-safe.
"""
function leading_term(polynomial; options...)
    first(leading_term([polynomial]; options...))
end

function leading_term(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:leading_term, options)
    result = leading_term0(polynomials, keywords)
    result
end

"""
    leading_ideal(polynomials; options...)

Returns generators of the ideal of the leading terms.

If the input is not a Groebner basis, computes a Groebner basis.

## Arguments

- `polynomials`: an array of polynomials. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl.

## Returns

- `basis`: the basis of the ideal of the leading terms.

## Possible Options

- `ordering`: Specifies the monomial ordering. Available monomial orderings are:
    - `InputOrdering()` for inferring the ordering from the given `polynomials`
      (default),
    - `Lex()` for lexicographic,
    - `DegLex()` for degree lexicographic,
    - `DegRevLex()` for degree reverse lexicographic,
    - `WeightedOrdering(weights)` for weighted ordering,
    - `ProductOrdering(args...)` for block ordering,
    - `MatrixOrdering(matrix)` for matrix ordering.
    For details and examples see the corresponding documentation page.

## Example

Using AbstractAlgebra.jl:

```jldoctest_
using Groebner, Nemo
R, (x, y) = QQ["x", "y"]
leading_ideal([x*y^2 + x, y*x^2 + y])
```

## Notes

- The function is thread-safe.
- For AbstractAlgebra.jl and Nemo.jl, the function is most efficient for
    polynomials over `GF(p)`, `Native.GF(p)`, and `QQ`.
"""
function leading_ideal(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:leading_ideal, options)
    result = leading_ideal0(polynomials, keywords)
    result
end

"""
    quotient_basis(polynomials; options...)

Returns a monomial basis of the quotient algebra of a zero-dimensional ideal.

If the input is not a Groebner basis, computes a Groebner basis.
If the input is not a zero-dimensional ideal, an error is raised.

## Arguments

- `polynomials`: an array of polynomials. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl.

## Returns

- `basis`: an array of monomials, a quotient basis.

## Possible Options

- `ordering`: Specifies the monomial ordering. Available monomial orderings are:
    - `InputOrdering()` for inferring the ordering from the given `polynomials`
      (default),
    - `Lex()` for lexicographic,
    - `DegLex()` for degree lexicographic,
    - `DegRevLex()` for degree reverse lexicographic,
    - `WeightedOrdering(weights)` for weighted ordering,
    - `ProductOrdering(args...)` for block ordering,
    - `MatrixOrdering(matrix)` for matrix ordering.
    For details and examples see the corresponding documentation page.

## Example

Using AbstractAlgebra.jl:

```jldoctest_
using Groebner, Nemo
R, (x, y) = QQ["x", "y"]
quotient_basis([x*y^2 + x, y*x^2 + y])
```

## Notes

- The function is thread-safe.
- For AbstractAlgebra.jl and Nemo.jl, the function is most efficient for
    polynomials over `GF(p)`, `Native.GF(p)`, and `QQ`.
"""
function quotient_basis(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:quotient_basis, options)
    result = quotient_basis0(polynomials, keywords)
    result
end

"""
    dimension(polynomials; options...)

Computes the (Krull) dimension of the ideal generated by `polynomials`.

If input is not a Groebner basis, computes a Groebner basis.

## Arguments

- `polynomials`: an array of polynomials. Supports polynomials from
    AbstractAlgebra.jl, Nemo.jl, and DynamicPolynomials.jl.

## Returns

- `dimension`: an integer, the dimension.

## Example

Using AbstractAlgebra.jl:

```jldoctest_
using Groebner, Nemo
R, (x, y) = QQ["x", "y"]
dimension([x*y^2 + x, y*x^2 + y])
```

## Notes

- The function is thread-safe.
- For AbstractAlgebra.jl and Nemo.jl, the function is most efficient for
    polynomials over `GF(p)`, `Native.GF(p)`, and `QQ`.
"""
function dimension(polynomials::AbstractVector; options...)
    keywords = KeywordArguments(:dimension, options)
    result = dimension0(polynomials, keywords)
    result
end
