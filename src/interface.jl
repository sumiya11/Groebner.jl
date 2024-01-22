# This file is a part of Groebner.jl. License is GNU GPL v2.

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
    - `:auto` for the automatic choice (default),
    - `:deterministic` for deterministic sparse linear algebra, 
    - `:randomized` for probabilistic sparse linear algebra.
- `threaded`: The use of multi-threading. Available options are: 
    - `:auto` for the automatic choice (default),
    - `:no` never use multi-threading, 
    - `:yes` allow the use of multi-threading.
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
  - `:no`, print nothing (default),
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

Or, say, in another monomial ordering:

```jldoctest
groebner([x*y^2 + x, y*x^2 + y], ordering=Lex(y, x))
```
"""
function groebner(polynomials::AbstractVector; options...)
    Base.require_one_based_indexing(polynomials)

    keywords = KeywordsHandler(:groebner, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    # NOTE: Type assertion is needed for type stability. This limits us to only
    # accept input arrays with concrete `eltype`
    result = _groebner0(polynomials, keywords)::typeof(polynomials)

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
end

"""
    groebner_learn(polynomials; options...)

Computes a Groebner basis of `polynomials` and emits the computation trace.

The trace can be used to speed up the computation of subsequent Groebner bases,
which should be specializations of the same ideal `groebner_learn` has been
applied to.

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

# Apply (same support, different coefficients)
flag, gb_2 = groebner_apply!(trace, [2x*y^2 + 3x, 4y*x^2 + 5y])

@assert flag
```

Using `groebner_learn` and `groebner_apply!` over different ground fields:

```jldoctest
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

Using `groebner_apply!` in batches (works only in `:degrevlex` at the moment):

```jldoctest
using Groebner, AbstractAlgebra
R, (x, y) = polynomial_ring(GF(2^31-1), ["x", "y"], ordering=:degrevlex)

# Learn
trace, gb_1 = groebner_learn([x*y^2 + x, y*x^2 + y])

# Create rings with some other moduli
R2, (x2, y2) = polynomial_ring(GF(2^30+3), ["x", "y"], ordering=:degrevlex)
R3, (x3, y3) = polynomial_ring(GF(2^27+29), ["x", "y"], ordering=:degrevlex)

# Two specializations of the same ideal
batch = ([2x2*y2^2 + 3x2, 4y2*x2^2 + 5y2], [4x3*y3^2 + 4x3, 5y3*x3^2 + 7y3])

# Apply for two sets of polynomials at once
flag, (gb_2, gb_3) = groebner_apply!(trace, batch)

@assert flag
@assert (gb_2, gb_3) == map(groebner, batch)
```

Perhaps, in a more involved example, we will compute Groebner bases of the
Katsura-9 system:

```jldoctest
using Groebner, AbstractAlgebra, BenchmarkTools

# Create the system
kat = Groebner.katsuran(9, k=ZZ, ordering=:degrevlex)

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

Observe the better amortized performance of the batched `groebner_apply!`.
"""
function groebner_learn(polynomials::AbstractVector; options...)
    Base.require_one_based_indexing(polynomials)

    keywords = KeywordsHandler(:groebner_learn, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    result = _groebner_learn0(polynomials, keywords)

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
end

"""
    groebner_apply!(trace, polynomials; options...)
    groebner_apply!(trace, batch::NTuple{N, Vector}; options...)

Computes a Groebner basis of `polynomials` using the given `trace`. The input
`polynomials` must be an array of polynomials over a finite field.

It is possible to input a tuple of `N` arrays to compute `N` Groebner bases
simultaneously, which could be more efficient than computing them separately.

This function is **not** thread-safe, since it mutates the `trace`.

See also `groebner_learn`.

## Possible Options

The `groebner_apply!` routine automatically inherits most of its parameters from
the given `trace`.

## Example

For examples, see the documentation of `groebner_learn`.
"""
function groebner_apply! end

# Specialization for a single input
function groebner_apply!(trace, polynomials::AbstractVector; options...)
    Base.require_one_based_indexing(polynomials)

    keywords = KeywordsHandler(:groebner_apply!, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    result =
        _groebner_apply0!(trace, polynomials, keywords)::Tuple{Bool, typeof(polynomials)}

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
end
# Specialization for a batch of inputs
function groebner_apply!(
    trace,
    batch::NTuple{N, T}; # deliberately not using ::Tuple{T, Vararg{T, Nminus1}}
    options...
) where {N, T <: AbstractVector}
    @assert N in (1, 2, 4, 8, 16) "The batch size must be one of the following: 1, 2, 4, 8, 16"
    all(Base.require_one_based_indexing, batch)

    keywords = KeywordsHandler(:groebner_apply!, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    result = _groebner_apply0!(trace, batch, keywords)::Tuple{Bool, typeof(batch)}

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
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
    Base.require_one_based_indexing(polynomials)

    keywords = KeywordsHandler(:isgroebner, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    result = _isgroebner0(polynomials, keywords)::Bool

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
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
    Base.require_one_based_indexing(basis)
    Base.require_one_based_indexing(to_be_reduced)

    keywords = KeywordsHandler(:normalform, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    result = _normalform0(basis, to_be_reduced, keywords)::typeof(to_be_reduced)

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
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
    Base.require_one_based_indexing(basis)

    keywords = KeywordsHandler(:kbase, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    result = _kbase0(basis, keywords)::typeof(basis)

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
end

"""
    fglm(basis, ordering_from, ordering_to; options...)

Converts a Groebner basis from one monomial ordering to another.
"""
function fglm(
    basis::AbstractVector,
    ordering_from::AbstractMonomialOrdering,
    ordering_to::AbstractMonomialOrdering;
    options...
)
    Base.require_one_based_indexing(basis)

    keywords = KeywordsHandler(:fglm, options)

    logging_setup(keywords)
    statistics_setup(keywords)

    result = _fglm0(basis, ordering_from, ordering_to, keywords)::typeof(basis)

    performance_counters_print(keywords)
    statistics_print(keywords)

    result
end
