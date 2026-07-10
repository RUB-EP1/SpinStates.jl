# SpinStates.jl

SpinStates.jl tracks the spin state of a particle through a decay program in either
the **helicity** or **canonical** basis, reporting how its spin-projection
coefficients evolve — phase and 4π spinor branch included — under any sequence of
frame transformations. It builds on
[InstructionalDecayTrees.jl](https://github.com/RUB-EP1/InstructionalDecayTrees.jl)
for the kinematics and
[PartialWaveFunctions.jl](https://github.com/JuliaHEP/PartialWaveFunctions.jl) for
the Wigner ``D``-matrix.

A state is ``|\text{FourVector}, 2s, \text{coeffs}\rangle``: the carrier
four-vector fixes the momentum geometry ``(\phi,\theta,\xi)``, and
``\text{coeffs}\in\mathbb{C}^{2s+1}`` carry the spin content plus the spinor phase.
Under a transform with accumulated ``SU(2)`` matrix ``U`` the coefficients rotate by
the **Wigner rotation**

```math
w = u_\text{basis}(p')^{-1}\, U\, u_\text{basis}(p),
\qquad \text{coeffs} \leftarrow D^{s}(w)\,\text{coeffs},
```

with prep ``u_H = R(\phi,\theta)\,B_z`` (helicity) or
``u_C = R\,B_z\,R^{-1}`` (canonical). ``w`` is built entirely in ``SU(2)`` — never
decoded from the ``SO(3)`` block, which would lose the ``\pm1`` branch.

## Conventions used in the tutorials

Every tutorial uses the same tiny toolkit. Rotations and boosts are the ``SU(2)`` /
``SL(2,\mathbb C)`` matrices

```math
U_\text{rot}(\hat n,\theta) = e^{-i\theta\,\hat n\cdot\vec\sigma/2},
\qquad
U_\text{boost}(\hat n,\xi) = e^{+\xi\,\hat n\cdot\vec\sigma/2},
```

are provided by SpinStates as [`Urot`](@ref) and [`Uboost`](@ref), and match the
`Rx`/`Ry`/`Rz`/`Bz` convention of
[FourVectors.jl](https://github.com/mmikhasenko/FourVectors.jl) (`Bz` is
parameterised by ``\gamma=\cosh\xi``):

```@example toolkit
using SpinStates, FourVectors, InstructionalDecayTrees, LinearAlgebra

U = Urot([0, 1, 0], 0.7)             # SU(2) rotation about ŷ
U ≈ InstructionalDecayTrees._su2_ry(0.7)
```

You rarely build these by hand: `Rx`, `Ry`, `Rz`, and `Bz` from FourVectors are
**extended to act on a [`SpinState`](@ref)** directly (Tutorials
[2](rotations.md) and [3](boosts.md)), transforming the carrier momentum and the
spin coefficients together.

## Quick example

```@example toolkit
p1 = FourVector(0.3, 0.2, 0.1; M = 0.2)
p2 = FourVector(-0.1, 0.4, -0.3; M = 0.4)
p3 = FourVector(-0.2, -0.6, 0.2; M = 0.3)
objs = (p1, p2, p3)

s = spin_state(Helicity(), p1, 1, 1 // 2)                    # |+1/2⟩ helicity
path = (ToHelicityFrame((1, 2, 3)), ToHelicityFrame((1, 2)))
(final_objs, results, s′) = track_spin(path, objs, 1, s)

s′
```

A [`SpinState`](@ref) renders as its ket expansion in the chosen basis. The tutorials
that follow demonstrate what the coefficients mean and how they move.
```@contents
Pages = ["spin_expectation.md", "rotations.md", "boosts.md", "wigner_rotation.md", "dirac_spinors.md"]
Depth = 1
```

## References

**Foundations** — the helicity/canonical formalism, Wigner rotations, and spinor
conventions used here are standard textbook material:

- A. D. Martin and T. D. Spearman, *Elementary Particle Theory*, North-Holland
  (1970) — [Inspire](https://inspirehep.net/literature/2104945).
- S. U. Chung, *Spin Formalisms* (updated version II), CERN Yellow Report
  CERN-71-8 — [CDS](https://cds.cern.ch/record/186421?ln=en).
- J. D. Richman, *An Experimenter's Guide to the Helicity Formalism*,
  CALT-68-1148 (1984) — [Inspire](https://inspirehep.net/literature/202987).
- É. Gourgoulhon, *Special Relativity in General Frames*, Springer (2013) —
  [Inspire](https://inspirehep.net/literature/1686144).

**Advanced applications** — how this basic bookkeeping is put to work in realistic
cascade decays:

- K. Habermann and M. Mikhasenko, *Wigner rotations for cascade reactions*,
  [Phys. Rev. D **111**, 056015 (2025)](https://inspirehep.net/literature/2827198)
  ([arXiv:2409.06913](https://arxiv.org/abs/2409.06913)).
- M. Mikhasenko *et al.*, *Dalitz-plot decomposition for three-body decays*,
  [Phys. Rev. D **101**, 034033 (2020)](https://inspirehep.net/literature/1758460)
  ([arXiv:1910.04566](https://arxiv.org/abs/1910.04566)).
