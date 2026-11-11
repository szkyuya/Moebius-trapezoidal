using MoebiusTrapezoidal
using Random
using Test

@testset "randomly shifted rule" begin
    @testset "half shift is the deterministic rule" begin
        x, w = moebius_nodes_shifted(81, 0.5)
        x0, w0 = moebius_nodes(81)
        @test x ≈ x0 atol = 1e-12
        @test w ≈ w0 rtol = 1e-12
    end

    @testset "zero shift skips infinity" begin
        x, w = moebius_nodes_shifted(81, 0.0)
        @test all(isfinite, x)
        @test sum(w) ≈ 1 rtol = 1e-6
    end

    @testset "every shift integrates smooth functions" begin
        for shift in (0.0, 0.1, 0.37, 0.9)
            x, w = moebius_nodes_shifted(243, shift)
            @test sum(w .* x .^ 2) ≈ 1 rtol = 1e-12
        end
    end

    @testset "unbiased with standard error" begin
        exact = sqrt(2 / π)
        r = moebius_trap_randomized(abs, 33; nshifts=64, rng=MersenneTwister(1))
        @test r.std_error > 0
        @test abs(r.value - exact) < 5 * r.std_error
        @test length(r.estimates) == 64
    end

    @testset "standard error decreases with n" begin
        e = [moebius_trap_randomized(abs, n; rng=MersenneTwister(2)).std_error for n in (16, 64, 256)]
        @test e[1] > e[2] > e[3]
    end

    @testset "reproducible with a seed" begin
        a = moebius_trap_randomized(abs, 17; rng=MersenneTwister(3))
        b = moebius_trap_randomized(abs, 17; rng=MersenneTwister(3))
        @test a.value == b.value && a.std_error == b.std_error
    end

    @testset "vector valued integrand" begin
        r = moebius_trap_randomized(x -> [x^2, abs(x)], 65; rng=MersenneTwister(4))
        @test length(r.value) == 2 && length(r.std_error) == 2
        @test r.value[1] ≈ 1 rtol = 1e-6
        @test abs(r.value[2] - sqrt(2 / π)) < 5 * r.std_error[2]
    end

    @testset "other weights, options and BigFloat" begin
        @test moebius_trap_randomized(x -> x^2, 243, :logistic; rng=MersenneTwister(5)).value ≈ π^2 / 3 rtol = 1e-6
        @test abs(moebius_trap_randomized(identity, 243; scale=2, loc=1, rng=MersenneTwister(6)).value) < 1e-8
        r = moebius_trap_randomized(x -> x^2, 729; T=BigFloat, rng=MersenneTwister(7))
        @test r.value isa BigFloat
        @test abs(r.value - 1) < big"1e-30"
    end

    @testset "argument checks" begin
        @test_throws ArgumentError moebius_trap_randomized(abs, 9; nshifts=1)
        @test_throws ArgumentError moebius_nodes_shifted(9, 1.0)
    end
end
