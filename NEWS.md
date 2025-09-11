## Version 0.10.0

- Interface change: keyword argument for `groebner` and friends is changed from `threaded=:yes/:no` to `tasks=N`. This allows to specify the number of tasks used in multi-threading. 
- Minor performance improvements.

## Version 0.9.2

- Switching to Documenter.jl for documentation.

## Version 0.9.0

- New exported functions: `quotient_basis`, `leading_ideal`, `dimension`.
- Renamed function `lead` to `leading_term`.  
