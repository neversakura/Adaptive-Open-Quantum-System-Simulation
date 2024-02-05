using OpenQuantumTools, OrdinaryDiffEq, DM, Printf, AVQD
include("argParser.jl")

nqbit = parsed_args["nqubit"]
lind_type = parsed_args["lind"] |> uppercase
gamma = parsed_args["gamma"]
gamma_p = parsed_args["gamma_p"]
gamma_m = parsed_args["gamma_m"]
pool_type = parsed_args["pool"]
rcut = parsed_args["rcut"]
rel = parsed_args["relative"]

inst_name = if lind_type == "Z"
    "test_benchmark_1/vectorized"
elseif lind_type == "PM"
    "test_benchmark_2/vectorized"
else
    throw(ArgumentError("Lindbladian type: $lind_type not supported."))
end

data_loc = "data"

entry = load_entry(data_loc, inst_name)

params = Dict(
    "prob" => "ASC",
    "type" => "vectorized",
    "nqubit" => nqbit,
    "tf" => 10,
    "dt" => 0.01,
    "lind" => lind_type,
    "gamma" => gamma,
    "gamma_p" => gamma_p,
    "gamma_m" => gamma_m,
    "rcut" => rcut,
    "pool" => pool_type,
    "adap" => rel ? "unrestricted" : "restricted"
)

println("Soving for vectorized solution with parameters:")
println(params)

tf = params["tf"]
dt = params["dt"]

Hd = -standard_driver(nqbit)
Hp = -alt_sec_chain(1, 0.5, 1, nqbit)

gamma = lind_type == "Z" ? gamma : [gamma_p, gamma_m]
Llist = if lind_type == "Z"
    build_z_ops(nqbit)
elseif lind_type == "PM"
    build_pm_ops(nqbit)
else
    throw(ArgumentError("Lindbladian type not supported."))
end

H = VectorizedEffectiveHamiltonian([(t) -> (1 - t / tf), (t) -> t / tf], [Hd, Hp], gamma, Llist)
u0 = ones(ComplexF64, 2^nqbit) |> normalize
ansatz = Ansatz(u0, relrcut=rcut, vectorize=true, pool=pool_type)
res = solve_avq(H, ansatz, [0, tf], dt)


save(entry, params, "data.jld2", "t", res.t, "state", res.u)
save(entry, params, "ansatz.jld2", "t", res.t, "theta", res.θ, "ansatz", res.A, "norm", res.norm)