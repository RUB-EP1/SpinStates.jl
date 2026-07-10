# Spin orientation ⟨S⟩

```@setup t1
using SpinStates, InstructionalDecayTrees, FourVectors, LinearAlgebra
σx = ComplexF64[0 1; 1 0]; σy = ComplexF64[0 -im; im 0]; σz = ComplexF64[1 0; 0 -1]; I2 = ComplexF64[1 0; 0 1]
ndotσ(n) = n[1]*σx + n[2]*σy + n[3]*σz
Urot(n, θ) = cos(θ/2)*I2 - im*sin(θ/2)*ndotσ(n)
Uboost(n, ξ) = cosh(ξ/2)*I2 + sinh(ξ/2)*ndotσ(n)
boost_x(p, ξ) = p |> Ry(-π/2) |> Bz(cosh(ξ)) |> Ry(π/2)
```

A spin state's coefficients are amplitudes in its **rest frame**. The observable that
says *where the spin points* is the expectation of the spin operator
``\hat S = \tfrac12\vec\sigma`` (for spin ``\tfrac12``),

```math
\langle\hat S\rangle = \frac{\langle\psi|\hat S|\psi\rangle}{\langle\psi|\psi\rangle},
```

a 3-vector on the Bloch sphere. [`spin_expectation`](@ref) computes it (and
[`spin_operators`](@ref) returns the ``(2s{+}1)``-dim spin matrices for any spin).

## Where it points

An eigenstate points along its quantization axis with length ``s``; a superposition
tilts away from it.

```@example t1
p = FourVector(0.3, 0.2, 0.1; M = 0.5)

up   = spin_state(Helicity(), p, ComplexF64[1, 0])            # |+1/2⟩
plus = spin_state(Helicity(), p, ComplexF64[1, 1] / sqrt(2)) # (|+⟩+|−⟩)/√2

(spin_expectation(up), spin_expectation(plus))
```

The first points along ``+\hat z`` with ``\|\langle S\rangle\| = \tfrac12``; the
second lies in the ``xy``-plane. A general complex superposition points anywhere on
the sphere:

```@example t1
c = normalize(ComplexF64[0.8, 0.3 + 0.5im])
v = spin_expectation(spin_state(Helicity(), p, c))
(v, norm(v))                # still length 1/2 — a pure state is fully polarized
```

## How it changes

Evolving the state moves ``\langle S\rangle`` rigidly: the coefficients rotate by the
Wigner rotation ``D^s(w)``, so ``\langle S\rangle`` rotates by the corresponding
``SO(3)`` rotation. Rotate the whole configuration about ``\hat y`` by ``\beta`` and
watch a canonical spin follow:

```@example t1
β = 0.7
p0 = FourVector(0.4, -0.2, 0.3; M = 0.4)
p1 = p0 |> Ry(β)                                 # momentum rotates
U  = Urot([0, 1, 0], β)                          # matching SU(2)

s0 = spin_state(Canonical(), p0, c)
s1 = evolve(s0, U, p1)

before = spin_expectation(s0)
after  = spin_expectation(s1)
Ry3 = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)] # SO(3) rotation about ŷ
after ≈ Ry3 * before
```

So ``\langle S\rangle`` is carried by the same rotation that moves the momentum — the
mean spin orientation transforms as an honest 3-vector. The next tutorials look at
*which* rotation the coefficients see, and how it differs between the two bases.
