# Connection to Dirac spinors

```@setup t5
using SpinStates, InstructionalDecayTrees, FourVectors, LinearAlgebra
σx = ComplexF64[0 1; 1 0]; σy = ComplexF64[0 -im; im 0]; σz = ComplexF64[1 0; 0 -1]; I2 = ComplexF64[1 0; 0 1]
ndotσ(n) = n[1]*σx + n[2]*σy + n[3]*σz
Urot(n, θ) = cos(θ/2)*I2 - im*sin(θ/2)*ndotσ(n)
Uboost(n, ξ) = cosh(ξ/2)*I2 + sinh(ξ/2)*ndotσ(n)
boost_x(p, ξ) = p |> Ry(-π/2) |> Bz(cosh(ξ)) |> Ry(π/2)
```

The two-component coefficients are the ``SL(2,\mathbb C)`` building blocks of a full
four-component **Dirac spinor**. In the Weyl (chiral) basis
(Peskin & Schroeder 3.50),

```math
u(p) = \begin{pmatrix}\sqrt{p\cdot\sigma}\,\xi\\[2pt]\sqrt{p\cdot\bar\sigma}\,\xi\end{pmatrix},
\qquad p\cdot\sigma = E\,\mathbb 1 - \vec p\cdot\vec\sigma,\quad
p\cdot\bar\sigma = E\,\mathbb 1 + \vec p\cdot\vec\sigma,
```

where ``\xi`` is a rest-frame 2-spinor — precisely the **canonical** coefficients
(spin quantized along the lab ``\hat z``).

```@example t5
pσ(p)  = p.E * I2 - (p.px * σx + p.py * σy + p.pz * σz)
pσb(p) = p.E * I2 + (p.px * σx + p.py * σy + p.pz * σz)

m  = 0.5
p  = FourVector(0.3, 0.2, 0.4; M = m)
ξ  = normalize(ComplexF64[0.6, 0.8])          # canonical 2-spinor

u = vcat(sqrt(Hermitian(pσ(p))) * ξ, sqrt(Hermitian(pσb(p))) * ξ)
round.(u; digits = 4)
```

## The lower block is the canonical preparation

``\sqrt{p\cdot\bar\sigma}`` is the Hermitian boost onto ``p`` — which is exactly
``\sqrt m`` times the canonical preparation ``u_C(p)=R\,B_z\,R^{-1}`` that SpinStates
uses:

```@example t5
sqrt(Hermitian(pσb(p))) ≈ sqrt(m) * SpinStates._prep_su2(Canonical(), p)
```

So `SpinState` coefficients *are* the ``\xi`` of a Dirac spinor, and `evolve` moves
them the same way the Dirac spinor moves. The construction also solves the Dirac
equation, ``(\not\!p - m)u = 0``:

```@example t5
Z = zeros(ComplexF64, 2, 2)
pslash = [Z pσ(p); pσb(p) Z]                  # γ^μ p_μ in the Weyl basis
norm((pslash - m * I) * u) < 1e-12
```

## Transformation is a matrix multiplication

Under a Lorentz transformation the Dirac spinor transforms by ``u \to S(\Lambda)u``,
block-diagonal in the Weyl basis: the right-handed block is our ``SU(2)/SL(2,\mathbb C)``
matrix ``U`` (the one `evolve` uses), the left-handed block is ``(U^\dagger)^{-1}``.
Boost along ``\hat x`` and compare:

```@example t5
U = Uboost([1, 0, 0], 0.5)                    # the transform's SL(2,C) matrix
S = [inv(U') Z; Z U]                          # S(Λ) = diag(S_L, S_R)
u_transformed = S * u

# rebuild the Dirac spinor from the *evolved* canonical coefficients
p′ = boost_x(p, 0.5)
s′ = evolve(spin_state(Canonical(), p, ξ), U, p′)
u_rebuilt = vcat(sqrt(Hermitian(pσ(p′))) * s′.coeffs, sqrt(Hermitian(pσb(p′))) * s′.coeffs)

u_rebuilt ≈ u_transformed
```

The two routes agree exactly: pushing the four-component spinor through ``S(\Lambda)``,
or evolving the two-component coefficients with SpinStates and re-assembling. The
Wigner rotation on the coefficients is the shadow of the ``SL(2,\mathbb C)`` matrix
acting on the full Dirac spinor.
