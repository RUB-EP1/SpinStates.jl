# API reference

```@docs
SpinStates
```

## States and bases

```@docs
AbstractSpinBasis
Helicity
Canonical
SpinState
spin_state
twos
projections
to_basis
```

## Transforms and evolution

`Rx`, `Ry`, `Rz`, and `Bz` from
[FourVectors.jl](https://github.com/mmikhasenko/FourVectors.jl) are extended to act
on a [`SpinState`](@ref) directly (e.g. `Ry(s, θ)` or `s |> Ry(θ)`), transforming the
carrier momentum and the spin coefficients together.

```@docs
Urot
Uboost
wigner_rotation
evolve
track_spin
```

## Observables

```@docs
spin_operators
spin_expectation
```
