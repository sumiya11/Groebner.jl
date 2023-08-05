## Groebner v0.4.3 Release notes

Versions of dependencies:
- Added depency on ExprTools.jl.

Changes in the interface **(may be breaking)**:

- The normalization of the output bases have changed. Just as before, the bases are normalized to have the leading coefficients equal to one. In this release, the leading coefficient is selected based on the ordering of the Groebner basis, instead of the ordering of the frontend polynomial implementation.

Other changes:

- Modified the internals of the F4 linear algebra.

- Introduced internal functions `performance_counters_enabled` and `@timed_block` for assessing performance.

## Groebner v0.4.2 Release notes

Fixed the bug with `groebner_apply!`, which crashed when called on a system over a finite field different from the one used in `groebner_learn`.

## Groebner v0.4.1 Release notes

Reinstate support for MultivariatePolynomials.jl version 0.4.

## Groebner v0.4.0 Release notes 

Minor bug fixes and improvements.

Versions of dependencies:
- Added depency on SIMD.jl.
- Updated version of MultivariatePolynomials.jl to 0.5.
- Updated version of AbstractAlgebra.jl to 0.31.

Changes in the interface:
- Keyword argument `rng` removed. Instead, the `seed` keyword argument is provided for setting the seed of internal rng.
- Changed keyword arguments `check`, `linalg`, `monoms`, and `loglevel`.
In particular, now it is possible to specify `monom = :sparse` to use a sparse monomial representation.
- Added functions `grobner_learn` and `groebner_apply!`.
