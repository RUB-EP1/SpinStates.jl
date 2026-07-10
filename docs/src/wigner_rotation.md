# 4. Wigner rotation and the axis–angle decomposition

```@setup t4
using SpinStates, FourVectors, LinearAlgebra
σx = ComplexF64[0 1; 1 0]; σy = ComplexF64[0 -im; im 0]; σz = ComplexF64[1 0; 0 -1]
boost_x(p, ξ) = p |> Ry(-π/2) |> Bz(cosh(ξ)) |> Ry(π/2)
```

Tutorials [2](rotations.md) and [3](boosts.md) showed that the coefficients of a spin
state are turned by *some* rotation when the carrier moves. This tutorial names that
rotation — the [Wigner rotation](https://en.wikipedia.org/wiki/Wigner_rotation) — and
ties it to the standard axis–angle picture.

## What `wigner_rotation` computes

For a state prepared at ``p_\text{from}`` and carried to ``p_\text{to}`` by a
transform with ``SU(2)`` matrix ``U``, the coefficients are multiplied by

```math
w \;=\; u_\text{basis}(p_\text{to})^{-1}\; U\; u_\text{basis}(p_\text{from}),
```

where ``u_\text{basis}(p)`` is the state-preparation matrix (the boost/rotation that
builds ``p`` from the rest frame). Read right to left, ``w`` takes the rest frame to
``p_\text{from}``, applies ``U`` to reach ``p_\text{to}``, then undoes the
preparation at ``p_\text{to}`` — a round trip **rest → rest**, so ``w`` is a pure
rotation even though ``U`` is a boost. That is the object [`wigner_rotation`](@ref)
returns, and `evolve` applies ``D^s(w)`` to the coefficients.

## Two boosts make a rotation

Start at rest, boost along ``\hat z`` (rapidity ``\xi_1``), then along ``\hat x``
(rapidity ``\xi_2``). Neither boost alone rotates a canonical spin, but their
composition does:

```@example t4
ξ1, ξ2 = 0.9, 0.7
rest = FourVector(0.0, 0.0, 0.0; M = 1.0)

p1 = rest |> Bz(cosh(ξ1))            # after boost 1 (along +z)
p2 = boost_x(p1, ξ2)                 # after boost 2 (along +x)

U2 = Uboost([1, 0, 0], ξ2)           # SL(2,C) of the second boost
w  = wigner_rotation(Canonical(), p1, p2, U2)
round.(w; digits = 4)
```

``w`` is unitary — the boosts cancelled, leaving a rotation.

## It is the unitary part of the boost product

There is a clean algebraic statement of the same fact: any
``SL(2,\mathbb C)`` element factors as a positive Hermitian matrix (a pure boost)
times a unitary (a rotation), ``U_2U_1 = H\,R``. That unitary ``R`` **is** the Wigner
rotation:

```@example t4
U1 = Uboost([0, 0, 1], ξ1)
Uprod = U2 * U1
H = sqrt(Hermitian(Uprod * Uprod'))   # pure-boost part
R = H \ Uprod                          # unitary (Wigner) part
w ≈ R
```

## Axis and angle

Writing ``w = \cos(\Omega/2)\,\mathbb 1 - i\sin(\Omega/2)\,\hat n\cdot\vec\sigma``,
the angle comes from the trace (``\operatorname{tr} w = 2\cos(\Omega/2)``) and the
axis from ``\operatorname{tr}(\sigma_k w)``:

```@example t4
Ω = 2 * acos(clamp(real(tr(w)) / 2, -1, 1))
n̂ = normalize(real.([im * tr(σ * w) for σ in (σx, σy, σz)]))
(Ω, round.(n̂; digits = 6))
```

The axis is ``\hat y``, perpendicular to the boost plane, as it must be. For two
**perpendicular** boosts the angle has the closed form
``\tan\Omega = \dfrac{\sinh\xi_1\,\sinh\xi_2}{\cosh\xi_1+\cosh\xi_2}``:

```@example t4
Ω ≈ atan(sinh(ξ1) * sinh(ξ2) / (cosh(ξ1) + cosh(ξ2)))
```

## The angle really rotates the spin

Finally, the same ``(\hat n,\Omega)`` — assembled into an ``SO(3)`` matrix by
[Rodrigues' formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) —
is precisely the rotation that moves the mean spin ``\langle S\rangle`` when the state
is boosted:

```@example t4
K  = [0 -n̂[3] n̂[2]; n̂[3] 0 -n̂[1]; -n̂[2] n̂[1] 0]      # [n̂]×
R3 = I + sin(Ω) * K + (1 - cos(Ω)) * K^2                # Rodrigues

s0 = spin_state(Canonical(), p1, normalize(ComplexF64[0.8, 0.3 + 0.5im]))
s1 = evolve(s0, U2, p2)
spin_expectation(s1) ≈ R3 * spin_expectation(s0)
```

So the abstract ``SU(2)`` Wigner rotation, the ``SL(2,\mathbb C)`` polar
decomposition, the closed-form angle, and the physical precession of
``\langle S\rangle`` are one object seen four ways. The little-group phase of the
helicity basis in [Tutorial 2](rotations.md) is the same construction with a diagonal
``w``.
