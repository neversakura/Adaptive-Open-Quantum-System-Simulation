function build_lind(op, gamma, nqbit)
    ops = [single_clause([op], [i], 1.0, nqbit) for i in 1:nqbit]
    [Lindblad(gamma, o) for o in ops]
end

build_z_lind(gamma, nqbit) = InteractionSet(build_lind("Z", gamma, nqbit)...)
build_z_ops(nqbit) = [TagOperator(single_clause(["Z"], [i], 1.0, nqbit), "Z"*string(i), nqbit) for i in 1:nqbit]

function build_pm_lind(γ₊, γ₋, nqbit)
    σp = (σx + 1im * σy)
    σm = (σx - 1im * σy)
    ops1 = [single_clause([σp], [i], 1.0, nqbit) for i in 1:nqbit]
    linds1 = [Lindblad(γ₊, o) for o in ops1]
    ops2 = [single_clause([σm], [i], 1.0, nqbit) for i in 1:nqbit]
    linds2 = [Lindblad(γ₋, o) for o in ops2]
    InteractionSet(linds1..., linds2...)
end

function build_pm_ops(nqbit)
    σp = (σx + 1im * σy)
    σm = (σx - 1im * σy)
    ops1 = [TagOperator(single_clause([σp], [i], 1.0, nqbit), "σ₊"*string(i), nqbit) for i in 1:nqbit]
    ops2 = [TagOperator(single_clause([σm], [i], 1.0, nqbit), "σ₋"*string(i), nqbit) for i in 1:nqbit]
    [ops1, ops2]
end

function get_energy(Hd, Hp, tf, t, sol::Vector{Matrix{T}}) where {T<:Number}
    H = Hamiltonian([(t) -> (1 - t / tf), (t) -> t / tf], [Hd, Hp], unit=:ħ, dimensionless_time=false)
    [real(tr(s * H(ti))) for (ti, s) in zip(t, sol)]
end

function get_energy(Hd, Hp, tf, t, sol::Vector{Vector{T}}) where {T}
    H = Hamiltonian([(t) -> (1 - t / tf), (t) -> t / tf], [Hd, Hp], unit=:ħ, dimensionless_time=false)
    num_samples = length(sol)
    num_points = length(t)
    avg_H = []
    std_H = []
    for i in 1:num_points
        tmp = zeros(num_samples)
        for j in 1:num_samples
            ψ = sol[j][i]
            tmp[j] = real(ψ' * H(t[i]) * ψ)
        end
        m = sum(tmp) / num_samples
        std = sqrt(sum((tmp .- m) .^ 2)) / num_samples
        push!(avg_H, m)
        push!(std_H, std)
    end
    avg_H, std_H
end

function get_rho(sol::Vector{Vector{T}}) where T
    nqbit = length(sol[1][1]) |> log2 |> Int
    num_samples = length(sol)
    ρ = []
    for i in eachindex(sol[1])
        temp = zeros(ComplexF64, 2^nqbit, 2^nqbit)
        for j in eachindex(sol)
            temp = temp + sol[j][i] * sol[j][i]'
        end
        push!(ρ, temp / num_samples)
    end
    ρ
end

function get_pop(Hd, Hp, tf, t, sol::Vector{Matrix{T}}; lvl=size(sol[1], 1), digits=8) where {T<:Number}
    H = Hamiltonian([(t) -> (1 - t / tf), (t) -> t / tf], [Hd, Hp], unit=:ħ, dimensionless_time=false)
    pop = []
    vp = nothing
    for i in eachindex(t)
        ρ = sol[i]
        w, v = eigen_decomp(H, t[i], lvl=lvl)
        if !(vp === nothing)
            dindx = find_degenerate(w, digits=digits)
            for idx in dindx
                v[:, idx] = (v[:, idx] / vp[:, idx])' * v[:, idx]
            end
        end
        vp = v
        push!(pop, real(diag(v'ρ * v)))
    end
    pop
end

function get_pop(Hd, Hp, tf, t, sol::Vector{Vector{T}}; lvl=length(sol[1][1]), digits=8) where T
    num_samples = length(sol)
    H = Hamiltonian([(t) -> (1 - t / tf), (t) -> t / tf], [Hd, Hp], unit=:ħ, dimensionless_time=false)
    pop = []
    vp = nothing
    mean = []
    std = []
    for i in eachindex(t)
        w, v = eigen_decomp(H, t[i], lvl=lvl)
        if !(vp === nothing)
            dindx = find_degenerate(w, digits=digits)
            for idx in dindx
                v[:, idx] = (v[:, idx] / vp[:, idx])' * v[:, idx]
            end
        end
        vp = v
        pop = []
        for s in sol
            st = s[i]
            push!(pop, abs2.(v' * st))
        end
        m = sum(pop)/num_samples
        push!(mean, m)
        push!(std, sqrt.(sum([(p .- m) .^ 2 for p in pop]))/num_samples)
    end
    mean, std
end

function trace_dist(ρ, σ)
    temp = ρ - σ 
    real(tr(sqrt(temp'*temp)))/2
end

to_mul_arr(data) = hcat(data...)'
ei_lab(nlvl) =hcat([latexstring(@sprintf("E_%d", i)) for i in 0:nlvl-1]...)