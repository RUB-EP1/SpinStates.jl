# SpinStates.jl

[![Test](https://github.com/RUB-EP1/SpinStates.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/RUB-EP1/SpinStates.jl/actions/workflows/Test.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://rub-ep1.github.io/SpinStates.jl/dev/)

SpinStates.jl tracks the spin state of a particle through a decay program in either
the **helicity** or **canonical** basis, reporting how its spin-projection
coefficients evolve — phase and 4π spinor branch included — under any sequence of
frame transformations. It builds on
[InstructionalDecayTrees.jl](https://github.com/RUB-EP1/InstructionalDecayTrees.jl)
for the kinematics/tracking and
[PartialWaveFunctions.jl](https://github.com/JuliaHEP/PartialWaveFunctions.jl) for
the Wigner D-matrix.

## Idea

A state is `|FourVector, 2s, coeffs⟩`: the carrier four-vector fixes the momentum
geometry `(ϕ, θ, ξ)`, and `coeffs ∈ ℂ^(2s+1)` carry the spin content plus the
spinor phase. Under a transform with accumulated `SU(2)` matrix `U`, the
coefficients rotate by the Wigner rotation

```
w = u_basis(p′)⁻¹ · U · u_basis(p),     coeffs ← Dˢ(w) · coeffs
```

with prep `u_H = R(ϕ,θ)·Bz` (helicity) or `u_C = R·Bz·R⁻¹` (canonical). `w` is
computed entirely in `SU(2)` — never decoded from the `SO(3)` block, which would
lose the ±1 branch — and lifted to spin `s` via `PartialWaveFunctions`.

## Example

```julia
using SpinStates, InstructionalDecayTrees, FourVectors

p1 = FourVector(0.3, 0.2, 0.1; M = 0.2)
p2 = FourVector(-0.1, 0.4, -0.3; M = 0.4)
p3 = FourVector(-0.2, -0.6, 0.2; M = 0.3)
objs = (p1, p2, p3)

# |+1/2⟩ helicity state on particle 1
s = spin_state(Helicity(), p1, 1, 1 // 2)

# evolve it through a boost chain; returns transformed objects, measurements, state
path = (ToHelicityFrame((1, 2, 3)), ToHelicityFrame((1, 2)))
(final_objs, results, s′) = track_spin(path, objs, 1, s)

s′.coeffs                     # helicity amplitudes in the final frame
to_basis(s′, Canonical())     # same state, canonical basis
```

Lower-level building blocks `wigner_rotation` and `evolve` are exported for custom
drivers.

## Notes

- Phase-exact for boost-based frames (`ToHelicityFrame`, `ToHelicityFrameParticle2`).
  Paths containing `PlaneAlign`/`ToGottfriedJacksonFrame` inherit the ±1 spinor
  branch of the IDT tracker's `U` for the rotation step; a small IDT-side change
  (an explicit `SU(2)` for `PlaneAlign`) makes those phase-exact too.
- Massive carriers only (the helicity/canonical prep needs a rest frame).
