# Randomly shifted Möbius-transformed trapezoidal rule with an error estimate.
#
# The integral of |x| with respect to the standard Gaussian weight, as in
# cot_trap.m, now with random shifts. The standard error is compared with the
# actual error. Run from the directory julia/MoebiusTrapezoidal with
#
#     julia --project=. examples/randomized.jl

using MoebiusTrapezoidal
using Printf
using Random

exact = sqrt(2 / π)
rng = MersenneTwister(2025)
nshifts = 16

@printf("%6s %14s %13s %11s\n", "n", "estimate", "actual error", "std error")
for m in 2:12
    n = 2^m
    r = moebius_trap_randomized(abs, n; nshifts, rng)
    @printf("%6d %14.10f %13.3e %11.3e\n", n, r.value, abs(r.value - exact), r.std_error)
end

println("\nEach row uses $nshifts random shifts, i.e. n * $nshifts function values.")
println("value ± 2 * std error is an approximate 95% confidence interval,")
println("so the actual error is mostly, but not always, below two standard errors.")
