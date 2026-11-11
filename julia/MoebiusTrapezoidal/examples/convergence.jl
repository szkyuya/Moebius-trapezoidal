# Convergence of the Möbius-transformed trapezoidal rule.
#
# This reproduces the example of the Matlab code cot_trap.m, the integral of
# |x|^p with respect to the standard Gaussian weight. Run from the directory
# julia/MoebiusTrapezoidal with
#
#     julia --project=. examples/convergence.jl

using MoebiusTrapezoidal
using Printf

p = 1  # parameter of the integrand |x|^p, odd p so that Γ((p+1)/2) is a factorial
exact = sqrt(2.0^p / π) * factorial((p - 1) ÷ 2)

ns = [2^m + 1 for m in 2:16]
errs = Float64[]
for n in ns
    value = moebius_trap(x -> abs(x)^p, n)
    push!(errs, abs(value - exact))
    @printf("n = %6d   estimate = %.10f   error = %.3e\n", n, value, errs[end])
end

# least squares slope of log(error) against log(n), as polyfit in cot_trap.m
lx, ly = log.(ns[3:end]), log.(errs[3:end])
mx, my = sum(lx) / length(lx), sum(ly) / length(ly)
slope = sum((lx .- mx) .* (ly .- my)) / sum((lx .- mx) .^ 2)
@printf("estimated rate of convergence: n^%.2f\n", slope)

# automatic refinement until a relative tolerance is met
r = moebius_quad(x -> abs(x)^p; rtol=1e-8)
@printf("moebius_quad: value = %.10f, error estimate = %.1e, n = %d\n", r.value, r.error, r.n)

# With Plots.jl installed, a convergence graph is obtained by
#     using Plots
#     plot(ns, errs; xscale=:log10, yscale=:log10, xlabel="Number of points",
#          ylabel="Absolute error", title="Convergence graph", legend=false)
