# Boosts

```@setup t3
using SpinStates, InstructionalDecayTrees, FourVectors, LinearAlgebra
σx = ComplexF64[0 1; 1 0]; σy = ComplexF64[0 -im; im 0]; σz = ComplexF64[1 0; 0 -1]; I2 = ComplexF64[1 0; 0 1]
ndotσ(n) = n[1]*σx + n[2]*σy + n[3]*σz
Urot(n, θ) = cos(θ/2)*I2 - im*sin(θ/2)*ndotσ(n)
Uboost(n, ξ) = cosh(ξ/2)*I2 + sinh(ξ/2)*ndotσ(n)
boost_x(p, ξ) = p |> Ry(-π/2) |> Bz(cosh(ξ)) |> Ry(π/2)
```

Under boosts the roles from the previous tutorial swap over.

## Helicity is invariant under a collinear boost

Helicity is the spin along the momentum; a boost **along that same direction** does
not change the momentum direction, so nothing happens to the coefficients — not even
a phase:

```@example t3
ξ  = acosh(1.8)                       # rapidity of the applied boost
pz = FourVector(0.0, 0.0, 0.5; M = 0.3)   # momentum along +z
c0 = normalize(ComplexF64[0.6, 0.8])

sh  = spin_state(Helicity(), pz, c0)
sh2 = evolve(sh, Uboost([0, 0, 1], ξ), pz |> Bz(cosh(ξ)))
sh2.coeffs ≈ c0
```

This is the property that makes helicity the convention of choice for sequential
boosts down a decay chain.

## A non-collinear boost rotates a canonical spin (Thomas–Wigner)

Boost a particle whose momentum is **not** along the boost axis and the canonical
state picks up a genuine rotation — the Thomas–Wigner rotation:

```@example t3
px  = FourVector(0.5, 0.0, 0.0; M = 0.3)   # momentum along +x
px2 = px |> Bz(cosh(ξ))                     # boosted along +z

w = wigner_rotation(Canonical(), px, px2, Uboost([0, 0, 1], ξ))
Ω = 2 * acos(clamp(real(tr(w)) / 2, -1, 1)) # rotation angle from tr w = 2cos(Ω/2)
```

The mean spin ``\langle S\rangle`` tilts by this angle in the ``xz``-plane:

```@example t3
s0 = spin_state(Canonical(), px, ComplexF64[1, 0])
s1 = evolve(s0, Uboost([0, 0, 1], ξ), px2)
(spin_expectation(s0), round.(spin_expectation(s1); digits = 4))
```

A boost is *not* a rotation, yet a sequence of them behaves as one — which is exactly
the [Wigner rotation](wigner_rotation.md) of the next tutorial. The helicity state,
meanwhile, would only rephase here: it is the boost **plane**, not the boost itself,
that the two bases see differently.
