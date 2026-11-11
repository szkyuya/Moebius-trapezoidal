# Root-mean-squared error of the two randomized rules.
#
# moebius_trap_randomized uses n points with a random shift. moebius_trap_rmt
# also draws the number of points at random from {⌊n/2⌋, …, n} (Definition 4.1
# of arXiv:2407.13650), which improves the worst-case RMSE by half an order.
#
# The comparison uses two integrands with the Gaussian weight:
#  * |x|, whose transformed Fourier coefficients decay regularly. Both rules
#    converge at the same rate n^-2, and moebius_trap_rmt has a larger
#    constant since it uses fewer points on average.
#  * cos(n θ(x)), where θ(x) is the inverse of x = -cot(θ/2), i.e. a frequency
#    equal to n after the transformation. This is the kind of integrand that
#    determines the worst case: the n point rule aliases it completely, so its
#    error does not decrease, while moebius_trap_rmt still converges.
#
# Run from the directory julia/MoebiusTrapezoidal with
#
#     julia --project=. examples/rmse_comparison.jl

using MoebiusTrapezoidal
using Printf
using Random

θ(x) = 2 * atan(1, -x)  # inverse of x = -cot(θ/2)
rmse(estimates, exact) = sqrt(sum(abs2, estimates .- exact) / length(estimates))

repetitions = 200
rng = MersenneTwister(2025)
ns = [32, 64, 128, 256, 512]

println("RMSE over $repetitions repetitions, Gaussian weight\n")
@printf("%5s | %11s %11s | %15s %13s\n", "n", "|x|, shift", "|x|, RMT", "cos(nθ), shift", "cos(nθ), RMT")
rows = Vector{Vector{Float64}}()
for n in ns
    exact_abs = sqrt(2 / π)
    f_osc(x) = cos(n * θ(x))
    exact_osc = moebius_trap(f_osc, 3^10)
    row = [
        rmse(moebius_trap_randomized(abs, n; nshifts=repetitions, rng).estimates, exact_abs),
        rmse(moebius_trap_rmt(abs, n; nsamples=repetitions, rng).estimates, exact_abs),
        rmse(moebius_trap_randomized(f_osc, n; nshifts=repetitions, rng).estimates, exact_osc),
        rmse(moebius_trap_rmt(f_osc, n; nsamples=repetitions, rng).estimates, exact_osc),
    ]
    push!(rows, row)
    @printf("%5d | %11.2e %11.2e | %15.2e %13.2e\n", n, row[1], row[2], row[3], row[4])
end

# least squares slope of log(RMSE) against log(n) for each column
lx = log.(ns)
mx = sum(lx) / length(lx)
slope(ly) = sum((lx .- mx) .* (ly .- sum(ly) / length(ly))) / sum((lx .- mx) .^ 2)
s = [slope(log.(getindex.(rows, j))) for j in 1:4]
@printf("%5s | %11.2f %11.2f | %15.2f %13.2f\n", "rate", s[1], s[2], s[3], s[4])
