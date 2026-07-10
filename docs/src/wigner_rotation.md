# Wigner rotation and the axis–angle decomposition

```@setup t4
using SpinStates, InstructionalDecayTrees, FourVectors, LinearAlgebra
σx = ComplexF64[0 1; 1 0]; σy = ComplexF64[0 -im; im 0]; σz = ComplexF64[1 0; 0 -1]; I2 = ComplexF64[1 0; 0 1]
ndotσ(n) = n[1]*σx + n[2]*σy + n[3]*σz
Urot(n, θ) = cos(θ/2)*I2 - im*sin(θ/2)*ndotσ(n)
Uboost(n, ξ) = cosh(ξ/2)*I2 + sinh(ξ/2)*ndotσ(n)
boost_x(p, ξ) = p |> Ry(-π/2) |> Bz(cosh(ξ)) |> Ry(π/2)
```

The composition of two non-collinear boosts is a boost **times a rotation** — the
[Wigner rotation](https://en.wikipedia.org/wiki/Wigner_rotation). This tutorial ties
the ``w`` that SpinStates applies to the coefficients to the standard axis–angle
picture.

Start at rest, boost along ``\hat z`` (rapidity ``\xi_1``), then along ``\hat x``
(rapidity ``\xi_2``):

```@example t4
ξ1, ξ2 = 0.9, 0.7
rest = FourVector(0.0, 0.0, 0.0; M = 1.0)

p1 = rest |> Bz(cosh(ξ1))    # after boost 1 (along +z)
p2 = boost_x(p1, ξ2)         # after boost 2 (along +x)

U1 = Uboost([0, 0, 1], ξ1)
U2 = Uboost([1, 0, 0], ξ2)
nothing # hide
```

## The Wigner rotation is the unitary part of the boost product

Any ``SL(2,\mathbb C)`` element factors as a positive Hermitian (pure boost) times a
unitary (rotation): ``U_2U_1 = H\,R``. That unitary ``R`` **is** the Wigner rotation,
and it is exactly what `wigner_rotation` returns for a canonical state:

```@example t4
Uprod = U2 * U1
H = sqrt(Hermitian(Uprod * Uprod'))   # pure-boost part
R = H \ Uprod                          # unitary (Wigner) part

w = wigner_rotation(Canonical(), p1, p2, U2)
w ≈ R
```

## Axis and angle

Writing ``w = \cos(\Omega/2)\,\mathbb 1 - i\sin(\Omega/2)\,\hat n\cdot\vec\sigma``, the
angle comes from the trace and the axis from ``\mathrm{tr}(\sigma_k w)``:

```@example t4
Ω = 2 * acos(clamp(real(tr(w)) / 2, -1, 1))
n̂ = normalize(real.([im * tr(σ * w) for σ in (σx, σy, σz)]))
(Ω, round.(n̂; digits = 6))
```

The axis is ``\hat y`` — perpendicular to the boost plane, as it must be. For two
**perpendicular** boosts the angle has the closed form
``\tan\Omega = \dfrac{\sinh\xi_1\,\sinh\xi_2}{\cosh\xi_1+\cosh\xi_2}``:

```@example t4
Ω ≈ atan(sinh(ξ1) * sinh(ξ2) / (cosh(ξ1) + cosh(ξ2)))
```

## The axis and angle really rotate the spin

As a sanity check, the same ``(\hat n,\Omega)`` — assembled into an ``SO(3)`` rotation
by [Rodrigues' formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula) —
is exactly the rotation that moves the mean spin ``\langle S\rangle`` when the state is
boosted:

```@example t4
K  = [0 -n̂[3] n̂[2]; n̂[3] 0 -n̂[1]; -n̂[2] n̂[1] 0]      # [n̂]×
R3 = I + sin(Ω) * K + (1 - cos(Ω)) * K^2                # Rodrigues

s0 = spin_state(Canonical(), p1, normalize(ComplexF64[0.8, 0.3 + 0.5im]))
s1 = evolve(s0, U2, p2)
spin_expectation(s1) ≈ R3 * spin_expectation(s0)
```

So the abstract ``SU(2)`` Wigner rotation, the ``SL(2,\mathbb C)`` polar decomposition,
the closed-form angle, and the physical precession of ``\langle S\rangle`` are all the
same object seen four ways.
