
"""
    groebner(polynomials; options...)

Computes a Groebner basis of `polynomials`.

The `polynomials` must be an array of polynomials from `AbstractAlgebra.jl`,
`Nemo.jl`, or `DynamicPolynomials.jl`. For `AbstractAlgebra.jl` and `Nemo.jl`,
the coefficients of polynomials must be either in `GF(p)` or in `QQ`.

This function is thread-safe. Groebner.jl may also use multi-threading
internally. See the option `threaded` below.

## Possible Options

The `groebner` routine takes the following optional arguments:
- `reduced`: If the returned basis must be autoreduced and unique (default is
  `true`).
- `ordering`: Specifies the monomial ordering. Available monomial orderings are: 
    - `InputOrdering()` for inferring the ordering from the given `polynomials`
      (default), 
    - `Lex(args...)` for lexicographic, 
    - `DegLex(args...)` for degree lexicographic, 
    - `DegRevLex(args...)` for degree reverse lexicographic, 
    - `WeightedOrdering(args...)` for weighted ordering, 
    - `ProductOrdering(args...)` for block ordering, 
    - `MatrixOrdering(args...)` for matrix ordering. For details and examples
  see the corresponding documentation page.
- `certify`: Certify the obtained basis. When this option is `false`, the
    algorithm is randomized, and the result is correct with high probability.
    When this option is `true`, the result is guaranteed to be correct in case
    the ideal is homogeneous (default is `false`). 
- `linalg`: Linear algebra backend. Available options are: 
    - `:auto`: for the automatic choice (default),
    - `:deterministic` for deterministic sparse linear algebra, 
    - `:randomized` for probabilistic sparse linear algebra.
- `threaded`: The use of multi-threading. Available options are: 
    - `:auto`: for the automatic choice (default),
    - `:no`: do not use multi-threading, 
    - `:yes`: use multi-threading.
    Additionally, you can set the environment variable `GROEBNER_NO_THREADED` to
    `1` to disable all multi-threading in Groebner.jl. In this case, the
    environment variable takes precedence over the `threaded` option.
- `monoms`: Monomial representation used in the computations. The algorithm
    tries to automatically choose the most suitable monomial representation.
    Otherwise, set `monoms` to one of the following: 
    - `:auto` for the automatic choice (default), 
    - `:dense` for classic dense exponent vectors,
    - `:packed` for packed representation, 
    - `:sparse` for sparse representation.
- `modular`: Modular computation algorithm. Only has effect when computing basis
    over the rational numbers. Possible options are:
    - `:auto` for the automatic choice (default),
    - `:classic_modular` for the classic multi-modular algorithm,
    - `:learn_and_apply` for the learn & apply algorithm.
- `seed`: The seed for randomization. Default value is `42`. Groebner uses
    `Random.Xoshiro` and `Random.MersenneTwister` for random number generation.
- `loglevel`: Logging level, an integer. Higher values mean less verbose.
    Default value is `0`, which means that only warnings are printed. Set
    `loglevel` to negative values, e.g., to `-3`, for debugging.
- `maxpairs`: The maximum number of critical pairs used at once in matrix
    construction in the F4 algorithm. By default, this is not specified. Tweak
    this option to try to lower the amount of RAM consumed.
- `homogenize`: Controls the use of homogenization in the algorithm. Possible
  options are:
    - `:yes`, always homogenize the input ideal,
    - `:no`, never homogenize the input ideal,
    - `:auto`, for the automatic choice (default).
- `statistics`: After the function exits, print some runtime statistics.
  Possible options are:
  - `:no`, do not print anything (default),
  - `:timings`, print the table with timings and allocations. Note that
    `Groebner.performance_counters_enabled()` must be set to `true` for this to
    have effect.

## Example

Using `DynamicPolynomials`:

```jldoctest
using Groebner, DynamicPolynomials
@polyvar x y
groebner([x*y^2 + x, y*x^2 + y])
```

Using `AbstractAlgebra`:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"]
groebner([x*y^2 + x, y*x^2 + y])
```

Or, in another monomial ordering:

```jldoctest
groebner([x*y^2 + x, y*x^2 + y], ordering=Lex(y, x))
```
"""
function groebner(polynomials::AbstractVector; options...)
    keywords = KeywordsHandler(:groebner, options)

    setup_logging(keywords)
    setup_statistics(keywords)

    # NOTE: Type assertion is needed for type stability. This limits us to only
    # accept input arrays with concrete `eltype`
    _groebner(polynomials, keywords)::typeof(polynomials)
end

"""
    groebner_learn(polynomials; options...)

Computes a Groebner basis of `polynomials` and emits the computation trace.

The trace can be then used to speed up the computation of subsequent Groebner
bases, which should be the specializations of the same ideal as the one
`groebner_learn` has been applied to.

The input `polynomials` must be an array of polynomials over a finite field.

This function is thread-safe.

See also `groebner_apply!`.

## Example

Using `groebner_learn` and `groebner_apply!` over the same ground field:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Apply (same ground field, different coefficients)
success, gb_2 = groebner_apply!(trace, [2x*y^2 + 3x, 4y*x^2 + 5y])

@assert success
```

Using `groebner_learn` and `groebner_apply!` over different ground fields:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Create a ring with a different modulo
_R, (_x, _y) = GF(2^30+3)["x", "y"]

# Apply (with a different modulo)
success, gb_2 = groebner_apply!(trace, [2_x*_y^2 + 3_x, 4_y*_x^2 + 5_y])

@assert success
@assert gb_2 == groebner([2_x*_y^2 + 3_x, 4_y*_x^2 + 5_y])
```
"""
function groebner_learn(polynomials::AbstractVector; options...)
    keywords = KeywordsHandler(:groebner_learn, options)

    setup_logging(keywords)
    setup_statistics(keywords)

    _groebner_learn(polynomials, keywords)
end

"""
    groebner_apply!(trace, polynomials; options...)

Computes a Groebner basis of `polynomials` using the given computation `trace`.

The input `polynomials` must be an array of polynomials over a finite field.

This function is *not* thread-safe.

See also `groebner_learn`.

## Example

Using `groebner_learn` and `groebner_apply!` over the same ground field:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Apply (same ground field, different coefficients)
success, gb_2 = groebner_apply!(trace, [2x*y^2 + 3x, 4y*x^2 + 5y])

@assert success
```

Using `groebner_learn` and `groebner_apply!` over different ground fields:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Create a ring with a different modulo
_R, (_x, _y) = GF(2^30+3)["x", "y"]

# Apply (with a different modulo)
success, gb_2 = groebner_apply!(trace, [2_x*_y^2 + 3_x, 4_y*_x^2 + 5_y])

@assert success
@assert gb_2 == groebner([2_x*_y^2 + 3_x, 4_y*_x^2 + 5_y])
```
"""
function groebner_apply! end

# Specialization for a single input
function groebner_apply!(trace, polynomials::AbstractVector; options...)
    keywords = KeywordsHandler(:groebner_apply!, options)

    setup_logging(keywords)
    setup_statistics(keywords)

    _groebner_apply!(trace, polynomials, keywords)::Tuple{Bool, typeof(polynomials)}
end

"""
    isgroebner(polynomials; options...)

Checks if `polynomials` forms a Groebner basis.

## Possible Options

The `isgroebner` routine takes the following optional arguments:
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
- `certify`: Use deterministic algorithm (default is `false`).
- `seed`: The seed for randomization. Default value is `42`. Groebner uses
    `Random.Xoshiro` and `Random.MersenneTwister` for random number generation.
- `loglevel`: Logging level, an integer. Higher values mean less verbose.
    Default value is `0`, which means that only warnings are printed. Set
    `loglevel` to negative values, e.g., `-1`, for debugging.
- `statistics`: After the function exits, print some runtime statistics.
    Possible options are:
    - `:no`, do not print anything (default),
    - `:timings`, print the table with timings and allocations. Note that
      `Groebner.performance_counters_enabled()` must be set to `true` for this to
      have effect.

## Example

Using `DynamicPolynomials`:

```jldoctest
using Groebner, DynamicPolynomials
@polyvar x y;
isgroebner([x*y^2 + x, y*x^2 + y])
```

Using `AbstractAlgebra`:

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = QQ["x", "y"]
isgroebner([x*y^2 + x, y*x^2 + y])
```
"""
function isgroebner(polynomials::AbstractVector; options...)
    keywords = KeywordsHandler(:isgroebner, options)

    setup_logging(keywords)
    setup_statistics(keywords)

    _isgroebner(polynomials, keywords)::Bool
end

"""
    normalform(basis, tobereduced; options...)

Computes the normal form of polynomials `tobereduced` w.r.t `basis`.
`tobereduced` can be either a single polynomial or an array of polynomials.

## Possible Options

The `normalform` routine takes the following optional arguments:
- `check`: Check if the given array `basis` forms a Groebner basis (default is `false`).
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
- `loglevel`: Logging level, an integer. Higher values mean less verbose.
  Default value is `0`, which means that only warnings are printed. Set
  `loglevel` to negative values, e.g., `-1`, for debugging.
- `statistics`: After the function exits, print some runtime statistics.
  Possible options are:
  - `:no`, do not print anything (default),
  - `:timings`, print the table with timings and allocations. Note that
    `Groebner.performance_counters_enabled()` must be set to `true` for this to
    have effect.

## Example

```jldoctest
using Groebner, DynamicPolynomials
@polyvar x y;
normalform([y^2 + x, x^2 + y], x^2 + y^2 + 1)
```
"""
function normalform(basis::AbstractVector, to_be_reduced::AbstractVector; options...)
    keywords = KeywordsHandler(:normalform, options)

    setup_logging(keywords)
    setup_statistics(keywords)

    _normalform(basis, to_be_reduced, keywords)::typeof(to_be_reduced)
end

normalform(basis::AbstractVector, to_be_reduced; options...) =
    first(normalform(basis, [to_be_reduced]; options...))

"""
    kbase(basis; options...)

Computes the basis of the quotient ideal of `basis` as a vector space.

## Possible Options

The `kbase` routine takes the following optional arguments:
- `check`: Check if the given array `basis` forms a Groebner basis (default is `false`).
- `loglevel`: Logging level, an integer. Higher values mean less verbose.
    Default value is `0`, which means that only warnings are printed. Set
    `loglevel` to negative values, e.g., `-1`, for debugging.
- `statistics`: After the function exits, print some runtime statistics.
    Possible options are:
    - `:no`, do not print anything (default),
    - `:timings`, print the table with timings and allocations. Note that
      `Groebner.performance_counters_enabled()` must be set to `true` for this to
      have effect.

## Example

```jldoctest
using Groebner, DynamicPolynomials
@polyvar x y;
kbase([y^2 + x, x^2 + y])
```
"""
function kbase(basis::AbstractVector; options...)
    keywords = KeywordsHandler(:kbase, options)

    setup_logging(keywords)
    setup_statistics(keywords)

    _kbase(basis, keywords)::typeof(basis)
end
