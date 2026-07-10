# 2. Rotations

```@setup t2
using SpinStates, FourVectors, LinearAlgebra
```

SpinStates extends the FourVectors transforms `Rx`, `Ry`, `Rz`, `Bz` so they act on
a [`SpinState`](@ref) as well as on a four-vector. Applying one moves the **carrier
momentum and the spin coefficients together**:

```@example t2
p0 = FourVector(0.4, -0.2, 0.3; M = 0.4)
s0 = spin_state(Canonical(), p0, normalize(ComplexF64[0.6, 0.8]))

s1 = Ry(s0, 0.6)     # rotate by 0.6 about ŷ; equivalently  s0 |> Ry(0.6)
s1
```

The two bases respond very differently to such a **pure rotation** ``R`` — this is
the cleanest way to tell them apart.

## Canonical: the spin follows, ``\text{coeffs}\to D(R)\,\text{coeffs}``

The canonical boost ``u_C(p)=R\,B_z\,R^{-1}`` is *covariant*
(``u_C(Rp)=R\,u_C(p)\,R^{-1}``), so the rotation the coefficients see is ``R`` itself
and they transform by ``D^s(R)`` — **independent of the momentum**. For spin ``\tfrac12``
that Wigner ``D``-matrix is just the ``SU(2)`` matrix [`Urot`](@ref):

```@example t2
β = 0.6
Ry(s0, β).coeffs ≈ Urot([0, 1, 0], β) * s0.coeffs
```

This is why canonical states are the natural home for ordinary angular-momentum
algebra: a rotation acts on the spin index and nothing else.

## Helicity: helicity is conserved (only a phase)

Helicity is the spin **along the momentum**. Under a rotation the momentum turns, and
the massive-particle little group carries the spin along with it — a rotation *about
the momentum*. Each helicity component therefore only picks up a phase; the
magnitudes (the probabilities of ``\lambda=+\tfrac12`` vs ``-\tfrac12``) are frozen:

```@example t2
sh = spin_state(Helicity(), p0, normalize(ComplexF64[0.6, 0.8]))
sh1 = Ry(sh, β)

abs.(sh1.coeffs) ≈ abs.(sh.coeffs)     # helicity is conserved
```

```@example t2
sh1.coeffs ./ sh.coeffs                # …the change is a pure phase per component
```

So under rotations: **canonical mixes** (by ``D(R)``), **helicity only rephases**.
Both preserve the norm and both agree with the rigid rotation of ``\langle S\rangle``
from [Tutorial 1](spin_expectation.md). The precise rotation the coefficients see —
and how to read off its axis and angle — is the subject of
[Tutorial 4](wigner_rotation.md).
