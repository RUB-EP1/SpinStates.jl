# 5. Connection to Dirac spinors

```@setup t5
using SpinStates, FourVectors, LinearAlgebra
σx = ComplexF64[0 1; 1 0]; σy = ComplexF64[0 -im; im 0]; σz = ComplexF64[1 0; 0 -1]; I2 = ComplexF64[1 0; 0 1]
boost_x(p, ξ) = p |> Ry(-π/2) |> Bz(cosh(ξ)) |> Ry(π/2)
```

The two-component coefficients are the ``SL(2,\mathbb C)`` building blocks of a full
four-component **Dirac spinor**. In the Weyl (chiral) basis
(Peskin & Schroeder Eq. 3.50), the ordering is left-handed on top, right-handed
below:

```math
u(p) = \begin{pmatrix}\psi_L\\ \psi_R\end{pmatrix}
     = \begin{pmatrix}\sqrt{p\cdot\sigma}\,\xi\\[2pt]\sqrt{p\cdot\bar\sigma}\,\xi\end{pmatrix},
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

## Which block is the canonical preparation?

Both blocks are canonical boosts of the *same* ``\xi`` — they differ only by the
**sign of the rapidity**, one for each chirality:

```math
\sqrt{p\cdot\sigma} = \sqrt{m}\,e^{-\xi\,\hat n\cdot\vec\sigma/2},
\qquad
\sqrt{p\cdot\bar\sigma} = \sqrt{m}\,e^{+\xi\,\hat n\cdot\vec\sigma/2}.
```

SpinStates uses the ``e^{+\xi\,\hat n\cdot\vec\sigma/2}`` convention (its `Uboost`
and the canonical preparation ``u_C``), which is the **right-handed** ``(0,\tfrac12)``
representation. In Peskin's ordering that sits in the **lower** block — hence it is
``\sqrt{p\cdot\bar\sigma}``, not ``\sqrt{p\cdot\sigma}``, that equals ``\sqrt m\,u_C``:

```@example t5
uC = SpinStates._prep_su2(Canonical(), p)      # canonical preparation e^{+ξ n·σ/2}
(sqrt(Hermitian(pσb(p))) ≈ sqrt(m) * uC,       # lower block  = √m · u_C
 sqrt(Hermitian(pσ(p)))  ≈ sqrt(m) * inv(uC))  # upper block  = √m · u_C⁻¹
```

(With the opposite Weyl ordering the two blocks — and the identification — simply
swap.) The construction also solves the Dirac equation, ``(\not\!p - m)u = 0``:

```@example t5
Z = zeros(ComplexF64, 2, 2)
pslash = [Z pσ(p); pσb(p) Z]                   # γ^μ p_μ in the Weyl basis
norm((pslash - m * I) * u) < 1e-12
```

## Transformation is a matrix multiplication

Under a Lorentz transformation the Dirac spinor transforms by ``u \to S(\Lambda)u``,
block-diagonal in the Weyl basis: the right-handed block is our
``SU(2)/SL(2,\mathbb C)`` matrix ``U`` (the one `evolve` uses), the left-handed block
is ``(U^\dagger)^{-1}``. Boost along ``\hat x`` and compare the two routes:

```@example t5
U = Uboost([1, 0, 0], 0.5)                     # the transform's SL(2,C) matrix
S = [inv(U') Z; Z U]                           # S(Λ) = diag(S_L, S_R)
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
