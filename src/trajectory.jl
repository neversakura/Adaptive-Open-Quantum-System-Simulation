# For distributed computing, setting --project won't set the package directory to the local folder
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using AVQD, OpenQuantumTools, DM

include("argParser.jl")

data_loc = "data"

nqbit = parsed_args["nqubit"]
rcut = parsed_args["rcut"]
rel = parsed_args["relative"]
lind = parsed_args["lind"] |> uppercase
gamma = parsed_args["gamma"]
gamma_p = parsed_args["gamma_p"]
gamma_m = parsed_args["gamma_m"]
pool_type = parsed_args["pool"]
postprocess = parsed_args["reduce"]
num_trajectory = parsed_args["num_t"]
no_state = parsed_args["no_state"]
save_every = parsed_args["save_every"]

@everywhere begin
    nqbit = $nqbit
    rcut = $rcut
    lind = $lind
    rel = $rel
    gamma = $gamma
    gamma_p = $gamma_p
    gamma_m = $gamma_m
    pool_type = $pool_type
    save_state = !($no_state)
    save_every = $save_every
    tf = 10
    dt = 0.01

    sp = nqbit > 5 ? true : false

    Hd = -standard_driver(nqbit, sp=sp)
    Hp = -alt_sec_chain(1, 0.5, 1, nqbit, sp=sp)

    gamma = lind == "Z" ? gamma : [gamma_p, gamma_m]
    Llist = if lind == "Z"
        build_z_ops(nqbit)
    elseif lind == "PM"
        build_pm_ops(nqbit)
    else
        throw(ArgumentError("Lindbladian type not supported."))
    end

    H = EffectiveHamiltonian([(t) -> (1 - t / tf), (t) -> t / tf], [Hd, Hp], gamma, Llist)
    u0 = ones(ComplexF64, 2^nqbit) |> normalize
    ansatz = rel ? Ansatz(u0, relrcut=rcut, pool=pool_type) : Ansatz(u0, rcut=rcut, pool=pool_type)

    struct Result
        t::Vector
        ψ::Vector
        H::Vector
        pop::Vector
        θ::Vector
        A::Vector
        jump_t::Vector
        jump_L::Vector
        status::Integer
    end

    function solve_avq_para(H, ansatz, tspan, dt, postprocess, save_state, save_every)
        res = solve_avq(H, ansatz, tspan, dt, save_state = save_state, save_everystep=save_every)
        tlist = res.t
        ψlist = res.u
        energy = []
        pop = []
        if postprocess && save_state
            vp = nothing
            for i in eachindex(tlist)
                Hm = AVQD.herm(H, tlist[i])
                ψ = ψlist[i]
                push!(energy, real(ψ' * Hm * ψ))

                w, v = eigen_decomp(H.Hh, tlist[i], lvl=10)
                if !(vp === nothing)
                    dindx = find_degenerate(w, digits=8)
                    for idx in dindx
                        v[:, idx] = (v[:, idx] / vp[:, idx])' * v[:, idx]
                    end
                end
                vp = v
                push!(pop, abs2.(v' * ψ))
            end
            # if postprocess is turned on, don't save the state
            ψlist = []
        end
        Result(tlist, ψlist, energy, pop, res.θ, res.A, res.jump_t, res.jump_L, res.status)
    end
end

params = Dict(
    "prob" => "ASC",
    "type" => "trajectory",
    "nqubit" => nqbit,
    "tf" => tf,
    "dt" => dt,
    "lind" => lind,
    "pool" => pool_type,
    "gamma" => gamma,
    "gamma_p" => gamma_p,
    "gamma_m" => gamma_m,
    "rcut" => rcut,
    "adap" => rel ? "unrestricted" : "restricted"
)
@printf "Solving %d AVQD trajectories with parameters:\n" num_trajectory
println(params)

inst_name = if lind == "Z"
    "test_benchmark_1/trajectory"
elseif lind == "PM"
    "test_benchmark_2/trajectory"
else
    throw(ArgumentError("Lindbladian type not supported."))
end

@printf "and instance name: %s.\n" inst_name

res = pmap((x) -> solve_avq_para(H, ansatz, [0, tf], dt, postprocess, save_state, save_every), 1:num_trajectory, batch_size=10)

entry = load_entry(data_loc, inst_name)

if save_state
    if postprocess
        save(entry, params, "plotData.jld2", "t", res[1].t, "H", [x.H for x in res], "pop", [x.pop for x in res], "status", [x.status for x in res])
    else
        save(entry, params, "data.jld2", "t", res[1].t, "state", [x.ψ for x in res], "status", [x.status for x in res])
    end
end
save(entry, params, "ansatz.jld2", "t", res[1].t, "theta", [x.θ for x in res], "ansatz", [x.A for x in res], "jump_t", [x.jump_t for x in res], "jump_L", [x.jump_L for x in res], "norm", [x.norm for x in res])