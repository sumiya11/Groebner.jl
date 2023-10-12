
"""
    groebner(polynomials; options...)

Computes a Groebner basis of `polynomials`.

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
    - `MatrixOrdering(args...)` for matrix ordering. 
  For details and examples see the corresponding documentation page.
- `certify`: Certify the obtained basis. When this option is `false`, the
    algorithm is randomized, and the result is correct with high probability.
    When this option is `true`, the result is guaranteed to be correct in case
    the ideal is homogeneous (default is `false`). 
- `linalg`: Linear algebra backend. Available options are: 
    - `:deterministic` for deterministic sparse linear algebra, 
    - `:randomized` for probabilistic sparse linear algebra (default is
      `:randomized`).
- `monoms`: Monomial representation used in the computations. The algorithm
    tries to automatically choose the most suitable monomial representation.
    Otherwise, set `monoms` to one of the following: 
    - `:auto` for the automatic choice (default), 
    - `:dense` for classic dense exponent vectors,
    - `:packed` for packed representation, 
    - `:sparse` for sparse representation.
- `seed`: The seed for randomization. Default value is `42`. Groebner uses
    `Random.Xoshiro` and `Random.MersenneTwister` for random number generation.
- `loglevel`: Logging level, an integer. Higher values mean less verbose.
    Default value is `0`, which means that only info and warnings are printed.
    Set `loglevel` to negative values, e.g., `-1`, for debugging.
- `maxpairs`: The maximum number of critical pairs used at once in matrix
    construction in the F4 algorithm. By default, this is not specified. Tweak
    this option to try to lower the amount of RAM consumed.
- `homogenize`: Controls the use of homogenization in the algorithm. Possible
  options are:
    - `:yes`, always homogenize the input ideal,
    - `:no`, never homogenize the input ideal,
    - `:auto`, for the automatic choice (default).

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
    # `KeywordsHandler` does several useful things on initialization:
    #   - checks that the keyword arguments are valid,
    #   - sets the global logger for this module,
    #   - refereshes performance counters.

    # NOTE: Type assertion is needed for type stability
    _groebner(polynomials, KeywordsHandler(:groebner, options))::typeof(polynomials)
end

"""
    groebner_learn(polynomials; options...)

Computes a Groebner basis of `polynomials` and emits the computation graph.

The graph can be used to speed up the computation of subsequent Groebner
bases, which should be the specializations of the same basis as the one
`groebner_learn` had been applied to.

*At the moment, only input over integers modulo a prime is supported.*

See also `groebner_apply!`.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
graph, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Apply
flag, gb_2 = groebner_apply!(graph, [2x*y^2 + 3x, 4y*x^2 + 5y])

@assert flag
```
"""
function groebner_learn(polynomials::AbstractVector; options...)
    _groebner_learn(polynomials, KeywordsHandler(:groebner_learn, options))
end

"""
    groebner_apply!(graph, polynomials; options...)

Computes a Groebner basis of `polynomials` using the given computation `graph`.

*At the moment, only input over integers modulo a prime is supported.*

See also `groebner_learn`.

## Example

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = GF(2^31-1)["x", "y"]

# Learn
graph, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Apply
flag, gb_2 = groebner_apply!(graph, [2x*y^2 + 3x, 4y*x^2 + 5y])

@assert flag
```
"""
function groebner_apply! end

# Specialization for a single input
function groebner_apply!(graph, polynomials::AbstractVector; options...)
    _groebner_apply!(
        graph,
        polynomials,
        KeywordsHandler(:groebner_apply!, options)
    )::Tuple{Bool, typeof(polynomials)}
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
    _isgroebner(polynomials, KeywordsHandler(:isgroebner, options))::Bool
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

## Example

```jldoctest
using Groebner, DynamicPolynomials
@polyvar x y;
normalform([y^2 + x, x^2 + y], x^2 + y^2 + 1)
```
"""
function normalform(basis::AbstractVector, tobereduced::AbstractVector; options...)
    _normalform(
        basis,
        tobereduced,
        KeywordsHandler(:normalform, options)
    )::typeof(tobereduced)
end

normalform(basis::AbstractVector, tobereduced; options...) =
    first(normalform(basis, [tobereduced]; options...))

"""
    kbase(basis; options...)

Computes the basis of the quotient ideal of `basis` as a vector space.

## Possible Options

The `kbase` routine takes the following optional arguments:
- `check`: Check if the given array `basis` forms a Groebner basis (default is `false`).
- `loglevel`: Logging level, an integer. Higher values mean less verbose.
    Default value is `0`, which means that only warnings are printed. Set
    `loglevel` to negative values, e.g., `-1`, for debugging.

## Example

```jldoctest
using Groebner, DynamicPolynomials
@polyvar x y;
kbase([y^2 + x, x^2 + y])
```
"""
function kbase(basis::AbstractVector; options...)
    _kbase(basis, KeywordsHandler(:kbase, options))::typeof(basis)
end
