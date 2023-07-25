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
