"""
    SpinStates

Track a particle's spin state through a decay program (built with
`InstructionalDecayTrees`) in either the **helicity** or **canonical** basis,
updating its projection coefficients — phase and 4π spinor branch included — under
any sequence of frame transformations.

A state is `|FourVector, 2s, coeffs⟩`: the carrier four-vector supplies the
momentum geometry `(ϕ, θ, ξ)`, while `coeffs ∈ ℂ^(2s+1)` hold the spin content and
the spinor phase. Under a transform with accumulated `SU(2)` matrix `U` the
coefficients rotate by the Wigner rotation `w = u_basis(p′)⁻¹·U·u_basis(p)`,
evaluated in `SU(2)` (never decoded from the `SO(3)` block, which loses the ±1
branch) and lifted to spin `s` via `PartialWaveFunctions`.

The conventions follow K. Habermann and M. Mikhasenko, *Wigner rotations for cascade
reactions*, Phys. Rev. D 111, 056015 (2025), arXiv:2409.06913; and M. Mikhasenko et
al., *Dalitz-plot decomposition for three-body decays*, Phys. Rev. D 101, 034033
(2020), arXiv:1910.04566.
"""
module SpinStates

using LinearAlgebra
using FourVectors: azimuthal_angle, polar_angle, boost_gamma
using PartialWaveFunctions: wignerD_doublearg
using InstructionalDecayTrees: InstructionalDecayTrees, init_tracked_state, apply_decay_instruction

const IDT = InstructionalDecayTrees

export AbstractSpinBasis, Helicity, Canonical
export SpinState, spin_state, projections, to_basis
export wigner_rotation, evolve, track_spin
export spin_operators, spin_expectation

"""
    AbstractSpinBasis

Supertype for the spin-quantization conventions: [`Helicity`](@ref) and
[`Canonical`](@ref). A basis selects the state-preparation element that maps a
rest-frame spin state onto the carrier momentum.
"""
abstract type AbstractSpinBasis end

"""
    Helicity()

Helicity basis: the state is prepared as `R(ϕ,θ)·Bz(ξ)|0,λ⟩` (boost along `+z`,
then rotate onto the momentum). The coefficient index is the helicity λ (spin
along the momentum).
"""
struct Helicity <: AbstractSpinBasis end

"""
    Canonical()

Canonical basis: the state is prepared as `R(ϕ,θ)·Bz(ξ)·R(ϕ,θ)⁻¹|0,m⟩`, the pure
(rotationless) boost onto the momentum. The coefficient index is the spin along
the fixed lab `z` axis.
"""
struct Canonical <: AbstractSpinBasis end

"""
    SpinState(basis, p, twos, coeffs)

Spin state of one particle carried by four-vector `p`, in `basis`
([`Helicity`](@ref) or [`Canonical`](@ref)). `twos = 2s`, `length(coeffs) = 2s+1`,
and `coeffs[k]` is projection `m = s - (k-1)` (index `1` is `+s`, `end` is `-s`;
see [`projections`](@ref)). Construct with [`spin_state`](@ref).
"""
struct SpinState{B <: AbstractSpinBasis, T <: Real, F}
    basis::B
    p::F
    twos::Int
    coeffs::Vector{Complex{T}}
end

"""
    spin_state(basis, p, coeffs)
    spin_state(basis, p, twos, m)

Build a [`SpinState`](@ref) on carrier four-vector `p`. The first form takes an
explicit coefficient vector (`twos = length(coeffs) - 1`); the second builds the
projection eigenstate `|m⟩` for spin `s = twos/2` (`m` integer or half-integer).
"""
function spin_state(basis::AbstractSpinBasis, p, coeffs::AbstractVector)
    cc = [complex(float(x)) for x in coeffs]
    T = real(eltype(cc))
    return SpinState{typeof(basis), T, typeof(p)}(basis, p, length(cc) - 1, cc)
end

function spin_state(basis::AbstractSpinBasis, p, twos::Integer, m::Real)
    twom = round(Int, 2m)
    (abs(twom) <= twos && iseven(twos - twom)) ||
        throw(ArgumentError("projection 2m=$twom incompatible with twos=$twos"))
    c = zeros(ComplexF64, twos + 1)
    c[(twos - twom) ÷ 2 + 1] = one(ComplexF64)
    return spin_state(basis, p, c)
end

"""
    projections(twos) -> Tuple

Projection values `m` (as `Rational{Int}`) aligned with the coefficient order of a
spin-`twos/2` [`SpinState`](@ref): `(s, s-1, …, -s)`.
"""
projections(twos::Integer) = ntuple(i -> (twos - 2 * (i - 1)) // 2, twos + 1)

# --- representation internals (reuse the SU(2)/decode convention of IDT) ---

# 2x2 SL(2,C) preparation elements built from the IDT SU(2) primitives. `p`
# supplies (ϕ, θ, ξ). Helicity uses two angles (the third ZYZ angle commutes with
# Bz and lives in `coeffs`); canonical is the rotationless boost R·Bz·R⁻¹.
function _prep_su2(::Helicity, p)
    ϕ = azimuthal_angle(p)
    θ = polar_angle(p)
    ξ = acosh(boost_gamma(p))
    return IDT._su2_rz(ϕ) * IDT._su2_ry(θ) * IDT._su2_bz(ξ)
end

function _prep_su2(::Canonical, p)
    ϕ = azimuthal_angle(p)
    θ = polar_angle(p)
    ξ = acosh(boost_gamma(p))
    return IDT._su2_rz(ϕ) * IDT._su2_ry(θ) * IDT._su2_bz(ξ) * IDT._su2_ry(-θ) * IDT._su2_rz(-ϕ)
end

# ZYZ Euler angles (ϕ, θ, ψ) reconstructing a unitary SU(2) matrix `w` as
# su2_rz(ϕ)·su2_ry(θ)·su2_rz(ψ). Angles are NOT wrapped: for half-integer spin the
# `e^{-imϕ}`/`e^{-imψ}` phases need the full 4π range, so wrapping ϕ into (-π,π]
# (as IDT's SO(3)-oriented decoder does) would flip the sign of Dˢ and break the
# homomorphism. Decoding from the SU(2) matrix (not the SO(3) block) keeps the
# spinor ±1 branch.
function _zyz_su2(w::AbstractMatrix)
    a, c = w[1, 1], w[2, 1]
    θ = 2 * atan(abs(c), abs(a))
    atol = 10 * eps(real(float(eltype(w))))
    if abs(c) < atol            # θ ≈ 0: only ϕ+ψ fixed
        return (ϕ = zero(real(θ)), θ = θ, ψ = 2 * angle(w[2, 2]))
    elseif abs(a) < atol        # θ ≈ π: only ϕ-ψ fixed
        return (ϕ = 2 * angle(c), θ = θ, ψ = zero(real(θ)))
    end
    σ = 2 * angle(w[2, 2])      # ϕ + ψ (mod 4π)
    δ = 2 * angle(c)            # ϕ - ψ (mod 4π)
    return (ϕ = (σ + δ) / 2, θ = θ, ψ = (σ - δ) / 2)
end

# (2s+1)-dim rotation representation of a *unitary* SU(2) matrix `w`, with elements
# from PartialWaveFunctions. Rows/cols run m = +s … -s.
function _wignerD(twos::Integer, w::AbstractMatrix)
    a = _zyz_su2(w)
    cosβ = cos(a.θ)
    ms = twos:-2:-twos
    return ComplexF64[
        wignerD_doublearg(twos, m1, m2, a.ϕ, cosβ, a.ψ) for m1 in ms, m2 in ms
    ]
end

"""
    wigner_rotation(basis, p_from, p_to, U) -> Matrix

The `2×2` `SU(2)` Wigner rotation `u_basis(p_to)⁻¹ · U · u_basis(p_from)` induced on
a `basis` spin state whose carrier moves from `p_from` to `p_to` under a transform
with accumulated `SU(2)` matrix `U`. Unitary (rest→rest); the boost parts cancel.
"""
wigner_rotation(basis::AbstractSpinBasis, p_from, p_to, U::AbstractMatrix) =
    _prep_su2(basis, p_to) \ (U * _prep_su2(basis, p_from))

"""
    evolve(s::SpinState, U, p_new) -> SpinState

Evolve `s` to carrier momentum `p_new` under a transform with accumulated `SU(2)`
matrix `U`, applying `Dˢ(wigner_rotation(...))` to the coefficients. `U` may be a
single step or a full path's accumulated matrix (the update telescopes).
"""
function evolve(s::SpinState, U::AbstractMatrix, p_new)
    w = wigner_rotation(s.basis, s.p, p_new, U)
    return SpinState(s.basis, p_new, s.twos, _wignerD(s.twos, w) * s.coeffs)
end

"""
    to_basis(s::SpinState, newbasis) -> SpinState

Re-express the same physical state at the same carrier momentum in `newbasis`,
rotating the coefficients by `Dˢ(u_new(p)⁻¹ · u_old(p))`.
"""
function to_basis(s::SpinState, newbasis::AbstractSpinBasis)
    w = _prep_su2(newbasis, s.p) \ _prep_su2(s.basis, s.p)
    return SpinState(newbasis, s.p, s.twos, _wignerD(s.twos, w) * s.coeffs)
end

"""
    spin_operators(twos) -> (Sx, Sy, Sz)

The three `(2s+1)×(2s+1)` spin matrices for spin `s = twos/2`, in the coefficient
ordering of a [`SpinState`](@ref) (`m = +s … -s`). For `twos = 1` these are `σ/2`.
"""
function spin_operators(twos::Integer)
    ms = [(twos - 2 * (i - 1)) / 2 for i in 1:(twos + 1)]   # +s … -s
    s = twos / 2
    n = twos + 1
    Sp = zeros(ComplexF64, n, n)   # raising operator
    for j in 1:n
        i = j - 1                   # m_i = m_j + 1
        i >= 1 && (Sp[i, j] = sqrt(s * (s + 1) - ms[j] * (ms[j] + 1)))
    end
    Sm = collect(Sp')
    return ((Sp + Sm) / 2, (Sp - Sm) / (2im), Diagonal(ComplexF64.(ms)) |> Matrix)
end

"""
    spin_expectation(s::SpinState) -> Vector{Float64}

The spin expectation vector `⟨Ŝ⟩ = ⟨ψ|Ŝ|ψ⟩ / ⟨ψ|ψ⟩` (with `Ŝ = σ/2` for spin-½)
in the state's rest frame, computed from its coefficients. A pure state has
`‖⟨Ŝ⟩‖ = s`; the direction is the spin's mean orientation on the Bloch sphere.
"""
function spin_expectation(s::SpinState)
    Sx, Sy, Sz = spin_operators(s.twos)
    c = s.coeffs
    nrm = real(c' * c)
    return [real(c' * S * c) / nrm for S in (Sx, Sy, Sz)]
end

"""
    track_spin(path, objs, idx, s::SpinState) -> (final_objs, results, s′)

Run an `InstructionalDecayTrees` instruction `path` on `objs` and evolve the spin
state `s` of particle `idx` (whose carrier must be `objs[idx]`) through the
accumulated Lorentz transform. Returns the transformed objects, the measurement
`results`, and the evolved [`SpinState`](@ref).

Phase-exact for boost-based frames (`ToHelicityFrame`, `ToHelicityFrameParticle2`).
Paths containing `PlaneAlign`/`ToGottfriedJacksonFrame` inherit the ±1 spinor
branch of the IDT tracker's `U` for the rotation step.
"""
function track_spin(path, objs, idx::Integer, s::SpinState)
    (tracked, results) = apply_decay_instruction(path, init_tracked_state(objs))
    return (tracked.objs, results, evolve(s, tracked.tracker.U, tracked.objs[idx]))
end

# --- pretty printing -------------------------------------------------------

_basis_name(::Helicity) = "helicity"
_basis_name(::Canonical) = "canonical"
_basis_tag(::Helicity) = "h"
_basis_tag(::Canonical) = "c"

_halfint(twon) = iseven(twon) ? string(twon ÷ 2) : string(twon, "/2")
_halfint_tex(twon) = iseven(twon) ? string(twon ÷ 2) : "\\tfrac{$twon}{2}"

_m_text(twom) = twom == 0 ? "0" : string(twom < 0 ? "-" : "+", _halfint(abs(twom)))
_m_tex(twom) = twom == 0 ? "0" : string(twom < 0 ? "-" : "+", _halfint_tex(abs(twom)))

function _coeff_str(z; digits = 3)
    re, im = round(real(z); digits), round(imag(z); digits)
    im == 0 && return string(re)
    re == 0 && return string(im, "i")
    return string("(", re, im < 0 ? "" : "+", im, "i)")
end

# join terms with sign-aware " + " / " - "
function _join_terms(parts)
    isempty(parts) && return "0"
    out = parts[1]
    for p in parts[2:end]
        out *= startswith(p, "-") ? " - " * p[2:end] : " + " * p
    end
    return out
end

function _ket_expansion(s::SpinState; tex::Bool)
    parts = map(enumerate(twos_range(s.twos))) do (k, twom)
        c = _coeff_str(s.coeffs[k])
        if tex
            "$c\\,\\left|$(_halfint_tex(s.twos)),$(_m_tex(twom))\\right\\rangle"
        else
            "$c |$(_halfint(s.twos)),$(_m_text(twom))⟩"
        end
    end
    return _join_terms(parts)
end

twos_range(twos) = twos:-2:-twos

# compact, e.g. inside a tuple/array
Base.show(io::IO, s::SpinState) =
    print(io, "SpinState(", _basis_name(s.basis), ", s=", _halfint(s.twos), ")")

function Base.show(io::IO, ::MIME"text/plain", s::SpinState)
    p = s.p
    r(x) = round(x; digits = 3)
    println(
        io, "SpinState · ", _basis_name(s.basis), " · s=", _halfint(s.twos),
        " · p=(", r(p.px), ", ", r(p.py), ", ", r(p.pz), "; ", r(p.E), ")"
    )
    return print(io, "  ", _ket_expansion(s; tex = false))
end

function Base.show(io::IO, ::MIME"text/latex", s::SpinState)
    return print(
        io, "\$\$|\\psi\\rangle_{\\mathrm{$(_basis_tag(s.basis))}} = ",
        _ket_expansion(s; tex = true), "\$\$"
    )
end

end # module SpinStates
