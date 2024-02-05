using OpenQuantumTools, OrdinaryDiffEq, DM, Printf, AVQD
include("argParser.jl")

nqbit = parsed_args["nqubit"]
lind_type = parsed_args["lind"] |> uppercase
gamma = parsed_args["gamma"]
gamma_p = parsed_args["gamma_p"]
gamma_m = parsed_args["gamma_m"]

inst_name = if lind_type == "Z"
    "test_benchmark_1/exact"
elseif lind_type == "PM"
    "test_benchmark_2/exact"
else
    throw(ArgumentError("Lindbladian type: $lind_type not supported."))
end

entry = load_entry("data", inst_name)

params = Dict(
    "prob" => "ASC",
    "type" => "exact",
    "nqubit" => nqbit,
    "tf" => 10,
    "dt" => 0.01,
    "lind" => lind_type,
    "gamma" => gamma,
    "gamma_p" => gamma_p,
    "gamma_m" => gamma_m
)

println("Soving for exact solution with parameters:")
println(params)

tf = params["tf"]
dt = params["dt"]

Hd = -standard_driver(nqbit)
Hp = -alt_sec_chain(1, 0.5, 1, nqbit)

linds = if lind_type == "Z"
    build_z_lind(gamma, nqbit)
elseif lind_type == "PM"
    build_pm_lind(gamma_p, gamma_m, nqbit)
end

H = Hamiltonian([(s) -> (1 - s), (s) -> s], [Hd, Hp], unit=:ħ)
u0 = ones(ComplexF64, 2^nqbit, 2^nqbit) / 2^nqbit
annealing = Annealing(H, u0, interactions=linds)

sol = solve_lindblad(annealing, tf, alg=Tsit5(), abstol=1e-6, reltol=1e-6, saveat=0:dt:tf)
save(entry, params, "data.jld2", "t", sol.t, "rho", sol.u)