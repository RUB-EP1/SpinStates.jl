# 3. Boosts

```@setup t3
using SpinStates, FourVectors, LinearAlgebra
```

Boosts act on a [`SpinState`](@ref) the same way — `Bz(s, γ)` (parameterised by
``\gamma=\cosh\xi``) moves the carrier momentum and the coefficients together. Under
boosts the roles from [Tutorial 2](rotations.md) swap over.

## Helicity is invariant under a collinear boost

Helicity is the spin along the momentum; a boost **along that same direction** does
not turn the momentum, so nothing happens to the coefficients — not even a phase:

```@example t3
pz = FourVector(0.0, 0.0, 0.5; M = 0.3)          # momentum along +z
sh = spin_state(Helicity(), pz, normalize(ComplexF64[0.6, 0.8]))

Bz(sh, 1.8).coeffs ≈ sh.coeffs                    # boost along +z: unchanged
```

This is the property that makes helicity the convention of choice for the sequential
boosts down a decay chain.

## A non-collinear boost rotates a canonical spin (Thomas–Wigner)

Boost a particle whose momentum is **not** along the boost axis and a canonical state
picks up a genuine rotation — the Thomas–Wigner rotation. Take a particle moving
along ``+\hat x`` and boost it along ``+\hat z``:

```@example t3
px = FourVector(0.5, 0.0, 0.0; M = 0.3)          # momentum along +x
sc = spin_state(Canonical(), px, ComplexF64[1, 0])  # pure |+1/2⟩ along lab ẑ

sc2 = Bz(sc, 1.8)
sc2                                               # the lower component is now populated
```

The mean spin ``\langle S\rangle`` tilts in the ``xz``-plane by the Wigner angle:

```@example t3
(spin_expectation(sc), round.(spin_expectation(sc2); digits = 4))
```

A boost is *not* a rotation, yet a canonical state responds to one as if rotated; a
helicity state in the same situation would only rephase. That "boost that acts like a
rotation" is exactly the [Wigner rotation](wigner_rotation.md) of the next tutorial,
where we extract its axis and angle.
