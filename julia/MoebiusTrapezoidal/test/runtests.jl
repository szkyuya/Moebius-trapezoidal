using MoebiusTrapezoidal
using Test

# least squares slope of log(err) against log(n)
function rate(ns, errs)
    lx, ly = log.(ns), log.(errs)
    mx, my = sum(lx) / length(lx), sum(ly) / length(ly)
    return sum((lx .- mx) .* (ly .- my)) / sum((lx .- mx) .^ 2)
end

@testset "MoebiusTrapezoidal" begin
    @testset "Gaussian moments" begin
        for (p, exact) in ((0, 1.0), (2, 1.0), (4, 3.0), (6, 15.0))
            @test moebius_trap(x -> x^p, 243) ≈ exact rtol = 1e-12
        end
        for p in (1, 3, 5)
            @test abs(moebius_trap(x -> x^p, 81)) < 1e-12
        end
    end

    @testset "other weights" begin
        @test moebius_trap(x -> 1.0, 243, logistic_weight) ≈ 1 rtol = 1e-9
        @test moebius_trap(x -> x^2, 729, :logistic) ≈ π^2 / 3 rtol = 1e-9
        @test moebius_trap(x -> exp(-x^2), 81, nothing) ≈ sqrt(π) rtol = 1e-12
        ρ(x) = exp(-(x - 1)^2 / 8) / sqrt(8π)  # density of N(1, 2^2)
        @test moebius_trap(identity, 243, ρ; scale=2, loc=1) ≈ 1 rtol = 1e-10
    end

    @testset "finite smoothness, rate n^-2 for |x|" begin
        exact = sqrt(2 / π)
        ns = [2^m + 1 for m in 4:11]
        errs = [abs(moebius_trap(abs, n) - exact) for n in ns]
        @test -2.2 < rate(ns, errs) < -1.8
    end

    @testset "scale and loc" begin
        for (scale, loc) in ((0.5, 0.0), (2.0, 1.0), (3.0, -1.0))
            @test moebius_trap(x -> x^2, 729; scale, loc) ≈ 1 rtol = 1e-10
        end
    end

    @testset "nodes and weights" begin
        x, w = moebius_nodes(81)
        @test issorted(x)
        @test all(>(0), w)
        @test sum(w) ≈ 1 rtol = 1e-8
        x2, w2 = moebius_nodes(81; prune=false)
        @test length(x2) == 81 && length(w2) == 81
    end

    @testset "vector valued integrand" begin
        @test moebius_trap(x -> [x^2, x^4], 243) ≈ [1, 3] rtol = 1e-12
    end

    @testset "BigFloat" begin
        v = moebius_trap(x -> x^2, 729; T=BigFloat)
        @test v isa BigFloat
        @test abs(v - 1) < big"1e-35"
    end

    @testset "automatic refinement" begin
        r = moebius_quad(x -> x^2; rtol=1e-13)
        @test r.converged
        @test r.value ≈ 1 rtol = 1e-13
        @test r.value ≈ moebius_trap(x -> x^2, r.n) rtol = 1e-14
        r2 = moebius_quad(abs; rtol=1e-8)
        @test r2.converged
        @test abs(r2.value - sqrt(2 / π)) < 1e-8
    end

    @testset "argument checks" begin
        @test_throws ArgumentError moebius_trap(abs, 0)
        @test_throws ArgumentError moebius_trap(abs, 9; scale=0)
        @test_throws ArgumentError moebius_trap(abs, 9, :unknown)
    end
end
