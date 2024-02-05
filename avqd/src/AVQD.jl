module AVQD

using Reexport
using Combinatorics, StatsBase

@reexport using LinearAlgebra, Printf, LaTeXStrings, OpenQuantumBase

include("utilities.jl")
include("pauliOps.jl")
include("effH.jl")
include("ansatz.jl")
include("solve.jl")

export EffectiveHamiltonian, Ansatz, VectorizedEffectiveHamiltonian
export pool_len, pool_tag, solve_avq, solve_vq, build_pool
export build_lind, build_z_lind, build_pm_lind, build_z_ops, build_pm_ops
export get_energy, get_rho, get_pop, trace_dist
export to_mul_arr, ei_lab

end # module AVQD
