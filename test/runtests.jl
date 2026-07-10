using SpinStates
using InstructionalDecayTrees
using FourVectors
using LinearAlgebra
using Test

const SS = SpinStates
const IDT = InstructionalDecayTrees

# Shared fixtures
const P1 = FourVector(0.3, 0.2, 0.1; M = 0.2)
const P2 = FourVector(-0.1, 0.4, -0.3; M = 0.4)
const P3 = FourVector(-0.2, -0.6, 0.2; M = 0.3)
const P4 = FourVector(0.25, -0.1, 0.35; M = 0.35)

# Accumulated SU(2) and transformed objects after running `path` on `objs`.
function run_path(path, objs)
    (t, _) = apply_decay_instruction(path, init_tracked_state(objs))
    return (t.tracker.U, t.objs)
end

@testset "SpinStates.jl" begin

    @testset "_wignerD representation" begin
        # The rotation representation is defined for unitary SU(2) inputs.
        w = IDT._su2_rz(0.3) * IDT._su2_ry(0.7) * IDT._su2_rz(0.5)
        w2 = IDT._su2_ry(1.1) * IDT._su2_rz(-0.4)

        @test SS._wignerD(1, w) ≈ w atol = 1.0e-12    # spin-1/2 rep is the matrix itself
        @test SS._wignerD(1, -w) ≈ -w atol = 1.0e-12  # opposite spinor branch preserved

        for twos in (1, 2, 3, 4)
            @test SS._wignerD(twos, Matrix{ComplexF64}(I, 2, 2)) ≈
                Matrix{ComplexF64}(I, twos + 1, twos + 1) atol = 1.0e-12
            @test SS._wignerD(twos, w * w2) ≈ SS._wignerD(twos, w) * SS._wignerD(twos, w2) atol =
                1.0e-11
            u = IDT._su2_rz(0.9) * IDT._su2_ry(2.1) * IDT._su2_rz(-1.3)   # unitary -> unitary
            D = SS._wignerD(twos, u)
            @test D * D' ≈ Matrix{ComplexF64}(I, twos + 1, twos + 1) atol = 1.0e-11
        end

        α = 0.5
        @test SS._wignerD(2, IDT._su2_rz(α)) ≈ Diagonal([cis(-α), 1, cis(α)]) atol = 1.0e-12

        β = 0.7   # spin-1 small-d against analytic Wigner d^1(β)
        c, s = cos(β), sin(β)
        d1 = [
            (1 + c) / 2 -s / sqrt(2) (1 - c) / 2
            s / sqrt(2) c -s / sqrt(2)
            (1 - c) / 2 s / sqrt(2) (1 + c) / 2
        ]
        @test SS._wignerD(2, IDT._su2_ry(β)) ≈ d1 atol = 1.0e-12
    end

    @testset "Urot/Uboost and direct-apply transforms" begin
        # matrix builders match the IDT SU(2) convention
        @test Urot([0, 1, 0], 0.7) ≈ IDT._su2_ry(0.7) atol = 1.0e-12
        @test Urot([0, 0, 1], 0.7) ≈ IDT._su2_rz(0.7) atol = 1.0e-12
        @test Uboost([0, 0, 1], 0.5) ≈ IDT._su2_bz(0.5) atol = 1.0e-12
        # Rz/Ry/Rx/Bz applied to a SpinState == evolve with the matching matrix
        c0 = normalize(ComplexF64[0.6, 0.8])
        s = spin_state(Canonical(), P1, c0)
        @test Ry(s, 0.6).coeffs ≈ evolve(s, Urot([0, 1, 0], 0.6), Ry(s.p, 0.6)).coeffs
        @test Rz(s, 0.6).coeffs ≈ evolve(s, Urot([0, 0, 1], 0.6), Rz(s.p, 0.6)).coeffs
        @test Rx(s, 0.6).coeffs ≈ evolve(s, Urot([1, 0, 0], 0.6), Rx(s.p, 0.6)).coeffs
        @test Bz(s, 1.8).coeffs ≈ evolve(s, Uboost([0, 0, 1], acosh(1.8)), Bz(s.p, 1.8)).coeffs
        # piping works and moves the carrier momentum
        @test (s |> Rz(0.6)).p.px ≈ Rz(s.p, 0.6).px
        # canonical spin follows a pure rotation by exactly D(R)
        @test Ry(s, 0.6).coeffs ≈ SS._wignerD(1, Urot([0, 1, 0], 0.6)) * c0 atol = 1.0e-11
    end

    @testset "spin_operators and spin_expectation" begin
        Sx, Sy, Sz = spin_operators(1)
        @test Sx ≈ [0 1; 1 0] / 2
        @test Sy ≈ [0 -im; im 0] / 2
        @test Sz ≈ [1 0; 0 -1] / 2
        # spin-1 satisfies the su(2) algebra [Sx,Sy] = i Sz
        Jx, Jy, Jz = spin_operators(2)
        @test Jx * Jy - Jy * Jx ≈ im * Jz atol = 1.0e-12
        # a pure state points with |⟨S⟩| = s; eigenstate along its quantization axis
        @test spin_expectation(spin_state(Helicity(), P1, 1, 1 // 2)) ≈ [0, 0, 0.5]
        @test spin_expectation(spin_state(Helicity(), P1, ComplexF64[1, 1] / sqrt(2))) ≈ [0.5, 0, 0] atol =
            1.0e-12
        @test norm(spin_expectation(spin_state(Canonical(), P2, 2, 1))) ≈ 1 atol = 1.0e-12
    end

    @testset "Spin state basics" begin
        s = spin_state(Helicity(), P1, 1, 1 // 2)
        @test twos(s) == 1 && s.coeffs == ComplexF64[1, 0]
        @test s isa SpinState{1}   # doubled spin is a type parameter
        @test spin_state(Canonical(), P1, 2, -1).coeffs == ComplexF64[0, 0, 1]
        @test projections(1) == (1 // 2, -1 // 2)
        @test projections(2) == (1 // 1, 0 // 1, -1 // 1)
        @test_throws ArgumentError spin_state(Helicity(), P1, 1, 1)     # 2m parity mismatch
        @test_throws ArgumentError spin_state(Helicity(), P1, 1, 3 // 2) # |m| > s
    end

    # -----------------------------------------------------------------------
    # Helicity basis: physical properties
    # -----------------------------------------------------------------------
    @testset "Helicity: collinear boost is invariant" begin
        # Carrier and boosted system both along +z: no Wigner rotation.
        q1 = FourVector(0.0, 0.0, 0.3; M = 0.2)
        q2 = FourVector(0.0, 0.0, 0.5; M = 0.4)
        c0 = normalize(ComplexF64[0.6, 0.8])
        (_, _, s2) = track_spin((ToHelicityFrame((1, 2)),), (q1, q2), 1, spin_state(Helicity(), q1, c0))
        @test s2.coeffs ≈ c0 atol = 1.0e-12
    end

    @testset "Helicity: pure rotation conserves helicity (diagonal Wigner rot.)" begin
        # Under a pure rotation the massive-helicity little group is Rz(γ), so the
        # Wigner rotation is diagonal and helicity magnitudes are preserved.
        (U, objs_f) = run_path((PlaneAlign(1, 2),), (P1, P2, P3))
        for twos in (1, 2)
            c0 = normalize(ComplexF64[1.0 + 0.2im, -0.3, 0.5im][1:(twos + 1)])
            w = wigner_rotation(Helicity(), P3, objs_f[3], U)
            s2 = evolve(spin_state(Helicity(), P3, c0), U, objs_f[3])
            @test w ≈ Diagonal(diag(w)) atol = 1.0e-12
            @test abs.(s2.coeffs) ≈ abs.(c0) atol = 1.0e-11
            @test norm(s2.coeffs) ≈ 1 atol = 1.0e-11
        end
    end

    # -----------------------------------------------------------------------
    # Canonical basis: physical properties
    # -----------------------------------------------------------------------
    @testset "Canonical: pure rotation acts as D(R), momentum-independent" begin
        # Defining property of the canonical boost: under a pure rotation the spin
        # follows by exactly D(R), regardless of the carrier momentum.
        (U, objs_f) = run_path((PlaneAlign(1, 2),), (P1, P2, P3))
        for twos in (1, 2)
            c0 = normalize(ComplexF64[0.5 + 0.1im, 0.8, -0.2][1:(twos + 1)])
            DR = SS._wignerD(twos, U)
            for carrier in (1, 3)
                s2 = evolve(spin_state(Canonical(), (P1, P2, P3)[carrier], c0), U, objs_f[carrier])
                @test s2.coeffs ≈ DR * c0 atol = 1.0e-11
            end
        end
    end

    @testset "Canonical: non-collinear boost induces mixing" begin
        q1 = FourVector(0.5, 0.0, 0.0; M = 0.2)   # along +x
        q2 = FourVector(0.0, 0.0, 0.6; M = 0.4)   # boost axis along +z (non-collinear)
        (_, _, s2) = track_spin((ToHelicityFrame((1, 2)),), (q1, q2), 1, spin_state(Canonical(), q1, 1, 1 // 2))
        @test abs(s2.coeffs[2]) > 1.0e-3           # lower component populated
        @test norm(s2.coeffs) ≈ 1 atol = 1.0e-12   # but still unitary
    end

    # -----------------------------------------------------------------------
    # Helicity <-> Canonical relations
    # -----------------------------------------------------------------------
    @testset "Relation: to_basis(H->C) = D(R(ϕ,θ))" begin
        for twos in (1, 2, 3)
            c0 = normalize(ComplexF64[0.7, -0.4 + 0.3im, 0.5, -0.2im][1:(twos + 1)])
            sh = spin_state(Helicity(), P2, c0)
            sc = to_basis(sh, Canonical())
            R = IDT._su2_rz(azimuthal_angle(P2)) * IDT._su2_ry(polar_angle(P2))
            @test sc.coeffs ≈ SS._wignerD(twos, R) * c0 atol = 1.0e-12
            @test to_basis(sc, Helicity()).coeffs ≈ c0 atol = 1.0e-12   # roundtrip
        end
    end

    @testset "Relation: bases coincide for momentum along +z" begin
        pz = FourVector(0.0, 0.0, 0.5; M = 0.3)   # θ = 0, R(ϕ,θ) = I
        c0 = normalize(ComplexF64[0.6, 0.8])
        @test to_basis(spin_state(Helicity(), pz, c0), Canonical()).coeffs ≈ c0 atol = 1.0e-12
    end

    @testset "Relation: evolve-in-each-basis then convert agree" begin
        (U, objs_f) = run_path((ToHelicityFrame((1, 2, 3)), ToHelicityFrame((1, 2))), (P1, P2, P3))
        c0 = normalize(ComplexF64[0.7, -0.4 + 0.3im])
        sh = spin_state(Helicity(), P1, c0)
        sc = to_basis(sh, Canonical())
        sh2 = evolve(sh, U, objs_f[1])
        sc2 = evolve(sc, U, objs_f[1])
        @test to_basis(sc2, Helicity()).coeffs ≈ sh2.coeffs atol = 1.0e-11
    end

    # -----------------------------------------------------------------------
    # Driver
    # -----------------------------------------------------------------------
    @testset "track_spin end-to-end is unitary (boost path & GJ)" begin
        c0 = normalize(ComplexF64[0.5 + 0.1im, 0.8])
        (_, _, s2) = track_spin(
            (ToHelicityFrame((1, 2, 3)), ToHelicityFrame((1, 2))),
            (P1, P2, P3), 1, spin_state(Helicity(), P1, c0),
        )
        @test norm(s2.coeffs) ≈ 1 atol = 1.0e-11

        (_, _, s3) = track_spin(
            (ToGottfriedJacksonFrame((2, 3, 4), 2, 4),),
            (P1, P2, P3, P4), 1, spin_state(Helicity(), P1, c0),
        )
        @test norm(s3.coeffs) ≈ 1 atol = 1.0e-10
    end
end
