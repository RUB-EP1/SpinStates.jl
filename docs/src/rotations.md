# Rotations

```@setup t2
using SpinStates, InstructionalDecayTrees, FourVectors, LinearAlgebra
σx = ComplexF64[0 1; 1 0]; σy = ComplexF64[0 -im; im 0]; σz = ComplexF64[1 0; 0 -1]; I2 = ComplexF64[1 0; 0 1]
ndotσ(n) = n[1]*σx + n[2]*σy + n[3]*σz
Urot(n, θ) = cos(θ/2)*I2 - im*sin(θ/2)*ndotσ(n)
Uboost(n, ξ) = cosh(ξ/2)*I2 + sinh(ξ/2)*ndotσ(n)
boost_x(p, ξ) = p |> Ry(-π/2) |> Bz(cosh(ξ)) |> Ry(π/2)
```

The two bases respond very differently to a **pure rotation** ``R`` — this is the
cleanest way to tell them apart.

Set up a carrier and a rotation about ``\hat y``:

```@example t2
β  = 0.6
p0 = FourVector(0.4, -0.2, 0.3; M = 0.4)
p1 = p0 |> Ry(β)             # rotated momentum
U  = Urot([0, 1, 0], β)      # matching SU(2)
c0 = normalize(ComplexF64[0.6, 0.8])
nothing # hide
```

## Canonical: the spin follows, ``\text{coeffs}\to D(R)\,\text{coeffs}``

The canonical boost ``u_C(p)=R\,B_z\,R^{-1}`` is *covariant*
(``u_C(Rp)=R\,u_C(p)\,R^{-1}``), so the Wigner rotation collapses to ``w=R`` itself.
The coefficients transform by ``D^s(R)`` — **independent of the momentum**:

```@example t2
sc = evolve(spin_state(Canonical(), p0, c0), U, p1)
sc.coeffs ≈ SpinStates._wignerD(1, U) * c0
```

`wigner_rotation` confirms the residual is exactly ``U``:

```@example t2
wigner_rotation(Canonical(), p0, p1, U) ≈ U
```

This is why canonical states are the natural home for ordinary angular-momentum
algebra: a rotation acts on the spin index and nothing else.

## Helicity: helicity is conserved (diagonal Wigner rotation)

For a helicity state the little group of a massive particle under rotations is a
rotation *about the momentum*, ``R_z(\gamma)``. The Wigner rotation is therefore
**diagonal** — each helicity component only picks up a phase, magnitudes are frozen:

```@example t2
wh = wigner_rotation(Helicity(), p0, p1, U)
round.(wh; digits = 4)
```

```@example t2
sh = evolve(spin_state(Helicity(), p0, c0), U, p1)
abs.(sh.coeffs) ≈ abs.(c0)      # |amplitudes| unchanged: helicity is conserved
```

So under rotations: **canonical mixes** (by ``D(R)``), **helicity only rephases**.
Both preserve the norm, and both agree with the rigid rotation of ``\langle S\rangle``
from [the previous tutorial](spin_expectation.md).
