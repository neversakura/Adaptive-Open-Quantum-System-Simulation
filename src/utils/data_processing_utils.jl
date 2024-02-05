# ============ load exact solution ================
function exact_param_trans(nqubit::Integer, lind::String)
    Dict(
        "prob" => "ASC",
        "type" => "exact",
        "nqubit" => nqubit,
        "tf" => 10,
        "dt" => 0.01,
        "lind" => lind,
        "gamma" => 0.01,
        "gamma_p" => 1e-2,
        "gamma_m" => 1e-3
    )
end

function load_exact_state(data_loc, nqubit::Integer, lind::String)
    params = exact_param_trans(nqubit, lind)
    inst_name = load_inst_name(params)
    entry = load_entry(data_loc, inst_name*"/exact")
    load(entry, params, "data.jld2", "t", "rho")
end

function load_exact_plot(data_loc, nqubit::Integer, lind::String)
    params = exact_param_trans(nqubit, lind)
    inst_name = load_inst_name(params)
    entry = load_entry(data_loc, inst_name*"/exact")
    try
        load(entry, params, "plot.jld2", "t", "H", "pop")
    catch
        @warn "Energy and population data not detected. Try loading the original data."
        (te, se) = load(entry, params, "data.jld2", "t", "rho")
        Hd = -standard_driver(params["nqubit"])
        Hp = -alt_sec_chain(1, 0.5, 1, params["nqubit"])
        exc_H = get_energy(Hd, Hp, params["tf"], te, se)
        pop_e = get_pop(Hd, Hp, params["tf"], te, se)
        pop_e = hcat(pop_e...)'
        save(entry, params, "plot.jld2", "t", te, "H", exc_H, "pop", pop_e)
        te, exc_H, pop_e
    end
end

function transfer_exact_plot(plot_loc, data_loc, nqubit::Integer, lind::String)
    # load the plotting data for the exact solution
    te, exc_H, pop_e = load_exact_plot(data_loc, nqubit, lind)
    # set the parameters for the plot data entry
    params = exact_param_trans(nqubit, lind)
    inst_name = load_inst_name(params)
    entry = load_entry(plot_loc, inst_name*"/exact")
    
    save(entry, params, "plot.jld2", "t", te, "H", exc_H, "pop", pop_e)
end
# ============ load vectorized solution ================
function vectorized_param_trans(nqubit::Integer, lind::String, adap::String, rcut::Number, pool::String)
    Dict(
        "prob" => "ASC",
        "type" => "vectorized",
        "nqubit" => nqubit,
        "tf" => 10,
        "dt" => 0.01,
        "lind" => lind,
        "gamma" => 0.01,
        "gamma_p" => 1e-2,
        "gamma_m" => 1e-3,
        "pool" => pool,
        "rcut" => rcut,
        "adap" => adap
    )
end

function load_vectorized_entry(data_loc, nqubit::Integer, lind::String, adap, rcut, pool)
    params = vectorized_param_trans(nqubit, lind, adap, rcut, pool)
    inst_name = load_inst_name(params)
    entry = load_entry(data_loc, inst_name*"/vectorized")
    entry, params
end

function load_vectorized_state(data_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="all2")
    entry, params = load_vectorized_entry(data_loc, nqubit, lind, adap, rcut, pool)
    load(entry, params, "data.jld2", "t", "state")
end

function load_vectorized_plot(data_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="all2")
    entry, params = load_vectorized_entry(data_loc, nqubit, lind, adap, rcut, pool)
    try
        load(entry, params, "plot.jld2", "t", "H", "pop")
    catch
        @warn "Energy and population data not detected. Try loading the original data."
            (t, s) = load(entry, params, "data.jld2", "t", "state")
            s = [i for i in s]
            Hd = -standard_driver(params["nqubit"])
            Hp = -alt_sec_chain(1, 0.5, 1, params["nqubit"])
            exc_H = get_energy(Hd, Hp, params["tf"], t, s)
            pop_e = get_pop(Hd, Hp, params["tf"], t, s)
            t, exc_H, hcat(pop_e...)'
    end
end

function load_vectorized_ansatz(data_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="all2")
    entry, params = load_vectorized_entry(data_loc, nqubit, lind, adap, rcut, pool)
    load(entry, params, "ansatz.jld2", "t", "theta", "ansatz", "norm")
end

function transfer_vectorized_plot(plot_loc, data_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="all2")
    entry, params = load_vectorized_entry(plot_loc, nqubit, lind, adap, rcut, pool)
    t, exc_H, pop_e = load_vectorized_plot(data_loc, nqubit, lind, adap=adap, rcut=rcut, pool=pool)
    save(entry, params, "plot.jld2", "t", t, "H", exc_H, "pop", pop_e)
end

function transfer_vectorized_data(des_loc, ori_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="all2")
    entry, params = load_vectorized_entry(des_loc, nqubit, lind, adap, rcut, pool)
    (t, s) = load_vectorized_state(ori_loc, nqubit, lind, adap=adap, rcut=rcut, pool=pool)
    save(entry, params, "data.jld2", "t", t, "state", s)
end

# ============ load trajectory solution ================
function trajectory_param_trans(nqubit::Integer, lind::String, adap::String, rcut::Number, pool::String)
    Dict(
        "prob" => "ASC",
        "type" => "trajectory",
        "nqubit" => nqubit,
        "tf" => 10,
        "dt" => 0.01,
        "lind" => lind,
        "gamma" => 0.01,
        "gamma_p" => 1e-2,
        "gamma_m" => 1e-3,
        "pool" => pool,
        "rcut" => rcut,
        "adap" => adap
    )
end

function load_traj_plot(data_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="neighbor")
    params = trajectory_param_trans(nqubit, lind, adap, rcut, pool)
    inst_name = load_inst_name(params)
    entry = load_entry(data_loc, inst_name * "/trajectory")
    try
        load(entry, params, "plot.jld2", "t", "avg_H", "avg_pop", "std_H", "std_pop")
    catch
        @warn "Energy and population data not detected. Try loading the original data."
        Hd = -standard_driver(params["nqubit"])
        Hp = -alt_sec_chain(1, 0.5, 1, params["nqubit"])
        (ta, avg_H, std_H, pop_a, std_a) = if check(entry, params, "data.jld2", "state")
            (ta, sa) = load(entry, params, "data.jld2", "t", "state");
            (avg_H, std_H, pop_a, std_a) = get_avg_H_pop(ta, sa, Hd, Hp, params["tf"])
            ta, avg_H, std_H, pop_a, std_a
        else
            (ta, Ha, popa) = load(entry, params, "plotData.jld2", "t", "H", "pop")
            (avg_H, std_H, pop_a, std_a) = get_avg_H_pop(ta, Ha, popa)
            ta, avg_H, std_H, pop_a, std_a
        end
        save(entry, params, "plot.jld2", "t", ta, "avg_H", avg_H, "std_H", std_H, "avg_pop", pop_a, "std_pop", std_a)
        ta, avg_H, pop_a, std_H, std_a
    end
end

function load_traj_ansatz(data_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="neighbor")
    params = trajectory_param_trans(nqubit, lind, adap, rcut, pool)
    inst_name = load_inst_name(params)
    entry = load_entry(data_loc, inst_name * "/trajectory")
    load(entry, params, "ansatz.jld2", "t", "theta", "ansatz", "jump_t", "jump_L")
end

function transfer_traj_plot(plot_loc, data_loc, nqubit::Integer, lind::String; adap="unrestricted", rcut=1e-4, pool="neighbor")
    ta, avg_H, pop_a, std_H, std_a = load_traj_plot(data_loc, nqubit, lind, adap=adap, rcut=rcut, pool=pool)
    
    params = trajectory_param_trans(nqubit, lind, adap, rcut, pool)
    inst_name = load_inst_name(params)
    entry = load_entry(plot_loc, inst_name * "/trajectory")
    
    save(entry, params, "plot.jld2", "t", ta, "avg_H", avg_H, "std_H", std_H, "avg_pop", pop_a, "std_pop", std_a)
end

# =================== utility functions ======================
function load_inst_name(params)
    if params["lind"]=="Z"
        "test_benchmark_1"
    elseif params["lind"]=="PM"
        "test_benchmark_2"
    else
        throw(ArgumentError("Lindbladian type not supported."))
    end
end

function get_energy_traj(Hd, Hp, tf, t, sol::Vector{Vector{T}}) where {T}
    H = Hamiltonian([(t) -> (1 - t / tf), (t) -> t / tf], [Hd, Hp], unit=:ħ, dimensionless_time=false)
    num_samples = length(sol)
    num_points = length(t)
    avg_H = []
    std_H = []
    res = []
    for j in 1:num_samples
        tmp = zeros(num_points)
        for i in 1:num_points
            ψ = sol[j][i]
            tmp[i] = real(ψ' * H(t[i]) * ψ)
        end
        push!(res, tmp)
    end
    res
end

function energy_comp(ext_loc, cmp_loc, params)
    (te, He, ) = load_exact(ext_loc, params) 
    cmp_type = params["type"]
    if cmp_type == "vectorized"
        (tv, Hv, ) = load_vectorized_plot(cmp_loc, params);
    else
        throw(ArgumentError("Type $cmp_type not supported."))
    end
    p = plot(tv, Hv, linewidth=2, label=L"\mathrm{AVQD}", xlabel=L"t", ylabel=L"\langle H(t) \rangle", 
    xguidefontsize=16, yguidefontsize=16, legendfontsize = 16,  tickfont = font(16, "Times"), legendposition=:right)
    plot!(te, He, linewidth=2, linestyle=:dash, label=L"\mathrm{Exact}")
    p
end

function get_avg_H_pop(ta, sa, Hd, Hp, tf)
    d = size(Hd, 1)
    lvl = d >= 8 ? 8 : d
    avg_H, std_H = get_energy(Hd, Hp, tf, ta, sa);
    pop_a, std_a = get_pop(Hd, Hp, tf, ta, sa);
    pop_a = (pop_a|>to_mul_arr)[:,1:lvl]
    std_a = (std_a|>to_mul_arr)[:,1:lvl];
    avg_H, std_H, pop_a, std_a
end

function get_avg_H_pop(ta, Ha, popa)
    d = size(Ha, 1)
    lvl = d >= 8 ? 8 : d
    num_samples = length(Ha)
    num_points = length(popa[1])
    avg_H = sum(Ha)/num_samples
    std_H = sqrt.(sum([(h - avg_H).^2 for h in Ha]))/length(Ha)
    
    pop_a = []
    std_a = []
    for i in 1:num_points
        res = zeros(length(popa[1][1]))
        for j in 1:num_samples
            res += popa[j][i]
        end
        mr = res/num_samples
        fill!(res, 0.0)
        for j in 1:num_samples
            res += (popa[j][i] - mr).^2
        end
        stdr = sqrt.(res)/num_samples
        push!(pop_a, mr)
        push!(std_a, stdr)
    end
    pop_a = (pop_a|>to_mul_arr)[:,1:lvl]
    std_a = (std_a|>to_mul_arr)[:,1:lvl]
    
    avg_H, std_H, pop_a, std_a
end;

# ============ loading circuit results from Siyuan ===============
function string_to_dict(s::AbstractString)
    # Remove leading and trailing braces
    s = strip(s, ['{', '}'])
    # Split on commas
    pairs = split(filter(x -> !isspace(x), s), ",")
    # Split key-value pairs on colons and convert to a dictionary
    Dict([split(pair, ":") |> (x) -> (strip(x[1], ['"','\''])|>string, parse(Float64, x[2])) for pair in pairs])
end

function load_circ_res(filename)
    circ_res = readlines(filename)
    circ_dict = Dict{String, Dict}()
    for i in 0:19
        loc_dict = Dict{String, Dict}()
        loc_dict["raw"] = circ_res[i*10+4] |> string_to_dict
        loc_dict["miti"] = circ_res[i*10+6] |> string_to_dict
        loc_dict["dd"] = circ_res[i*10+8] |> string_to_dict
        loc_dict["all"] = circ_res[i*10+10] |> string_to_dict
        circ_dict[string(i*50+50)] = loc_dict
    end
    circ_dict
end
