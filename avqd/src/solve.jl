"""
    function solve_avq(H::EffectiveHamiltonian, A::Ansatz, tspan, dt::Number; save_state::Bool, save_everystep::Bool)

Solves one open-system trajectory using adaptive variational algorithm.

# Arguments
- `H::EffectiveHamiltonian`: `EffectiveHamiltonian` object of the open system model.
- `A::Ansatz`: `Ansatz` object of the variatonal simulation.
- `tspan`: a tuple like object (t0, tf) defining the initial time (t0) and final time (tf).
- `dt::Number`: time stepsize for the simulation.
- `save_state::Bool`: whether to save the state of the quantum system, defaults to true.
- `save_everystep::Bool`: whether to save the result at every time step or only at jump points, defaults to false.
"""
function solve_avq(H::EffectiveHamiltonian, A::Ansatz, tspan, dt::Number;
    save_state::Bool=true, save_everystep::Bool=false)
    # store inital ref state for reinitializtion
    ref_init = A |> get_ref
    # reset current state
    update_state!(A)
    # initial time
    t = tspan[1]
    # results recorder
    t_list = []
    u_list = []
    θ_list = []
    A_list = []
    jump_L_list = []
    jump_t_list = []
    # break flag
    break_flag = false
    # jump recorder
    Γ = 0
    q = rand(Float64)
    if save_everystep
        push!(t_list, t)
        push!(θ_list, A.θ |> copy)
        push!(A_list, [i |> tag for i in A.A])
        save_state ? push!(u_list, A |> get_state) : nothing
    end
    while t + dt <= tspan[2]
        @debug "Solving for t:" t
        He = herm(H, t)
        Ha = antiherm(H)
        @debug "Current ansats:" [i |> tag for i in A.A]
        if exp(-Γ) > q
            try
                one_step!(A, He, Ha, dt)
            catch
                break_flag = true
                break
            end
            ψ = A |> get_state
            Γ = Γ + 2 * real(ψ' * Ha * ψ) * dt
            if save_everystep
                push!(t_list, t + dt)
                push!(θ_list, A.θ |> copy)
                push!(A_list, [i |> tag for i in A.A])
            end
        else
            @debug "Applying jump operator"
            ψ = A |> get_state
            w = Weights([real(ψ' * LdL * ψ) for LdL in get_LdL(H)])
            L = sample(get_L(H), w)
            ψₙ = normalize(L * ψ)
            # record jump time
            push!(jump_t_list, t)
            # record jump operator
            push!(jump_L_list, L |> tag)
            # record time step for the ansatz
            push!(t_list, t + dt)
            # set new reference state
            set_ref!(A, ψₙ)
            # reset Ansatz
            reset!(A)
            # record the new ansatz parameters after jump
            # record θ
            push!(θ_list, A.θ |> copy)
            # record ansatz operators
            push!(A_list, [i |> tag for i in A.A])
            # reset jump recorder
            Γ = 0
            q = rand(Float64)
        end
        save_state ? push!(u_list, A |> get_state) : nothing
        t += dt
    end
    if !break_flag
        # record final ansatz operators and θ
        if (isempty(jump_t_list) || !isapprox(t, jump_t_list[end])) && !save_everystep
            push!(θ_list, A.θ |> copy)
            push!(A_list, [i |> tag for i in A.A])
        end
    end
    # reset ansatz to the inital condition
    set_ref!(A, ref_init)
    reset!(A)
    # construct solution object
    AVQDSol(t_list, u_list, θ_list, A_list, jump_L_list, jump_t_list, status=!break_flag ? 0 : 1)
end

"""
    function solve_avq(H::VectorizedEffectiveHamiltonian, A::Ansatz, tspan, dt::Number)

Solves the vectorized Lindblad ME using the adaptive variational algorithm.
"""
function solve_avq(H::VectorizedEffectiveHamiltonian, A::Ansatz, tspan, dt::Number)
    # store inital ref state for reinitializtion
    ref_init = A |> get_ref
    nqbit = A.nqbit
    # reset current state
    update_state!(A)
    # initial time
    t = tspan[1]
    # results recorder
    t_list = []
    u_list = []
    θ_list = []
    A_list = []
    norm_list = []
    push!(t_list, t)
    push!(u_list, reshape(A |> get_state, 2^nqbit, 2^nqbit))
    push!(θ_list, A.θ |> copy)
    push!(A_list, [i |> tag for i in A.A])
    push!(norm_list, 1.0)
    Γ = 0
    t_idx = 2
    while t + dt <= tspan[2]
        @debug "Solving for t:" t
        He = herm(H, t)
        Ha = antiherm(H)
        @debug "Current ansats:" [i |> tag for i in A.A]

        one_step!(A, He, Ha, dt)
        ψ = A |> get_state
        Γ = Γ + 2 * real(ψ' * Ha * ψ) * dt
        ρ = reshape(ψ, 2^nqbit, 2^nqbit)
        # renormalize density matrix
        ρ = ρ / tr(ρ)

        t_idx += 1
        t += dt
        push!(t_list, t)
        push!(u_list, ρ)
        push!(θ_list, A.θ |> copy)
        push!(A_list, [i |> tag for i in A.A])
        # record the norm shrinking
        push!(norm_list, exp(-Γ))
    end
    # reset ansatz to the inital condition
    set_ref!(A, ref_init)
    reset!(A)
    AVQDSol(t_list, u_list, θ_list, A_list, [], [], norm=norm_list)
end

"""
    function solve_avq(H::EffectiveHamiltonian, A::Ansatz, tspan, dt::Number, save_state::Bool=true, save_everystep::Bool=false)

Solves one open-system trajectory using non-adaptive variational algorithm. This
function is for test purpose.
"""
function solve_vq(H::EffectiveHamiltonian, A::Ansatz, tspan, dt::Number;
    save_state::Bool=true, save_everystep::Bool=false)
    # store inital ref state for reinitializtion
    ref_init = A |> get_ref
    # set the ansazt to be every operator in the pool
    set_pool_to_ansatz!(A)
    # reset current state
    update_state!(A)
    # initial time
    t = tspan[1]
    # results recorder
    t_list = []
    u_list = []
    θ_list = []
    A_list = [i |> tag for i in A.A]
    jump_L_list = []
    jump_t_list = []
    # jump recorder
    Γ = 0
    q = rand(Float64)
    if save_everystep
        push!(t_list, t)
        push!(θ_list, A.θ |> copy)
        save_state ? push!(u_list, A |> get_state) : nothing
    end
    t_idx = 1
    while t + dt <= tspan[2]
        #@info "Solving for t:" t
        He = herm(H, t)
        Ha = antiherm(H)
        #@info "Current ansats:" [i|>tag for i in A.A]
        if exp(-Γ) > q
            one_step_non_adp!(A, He, Ha, dt)
            ψ = A |> get_state
            Γ = Γ + 2 * real(ψ' * Ha * ψ) * dt
            if save_everystep
                push!(t_list, d + dt)
                push!(θ_list, A.θ |> copy)
            end
        else
            #@info "Applying jump operator"
            ψ = A |> get_state
            w = Weights([real(ψ' * LdL * ψ) for LdL in get_LdL(H)])
            L = sample(get_L(H), w)
            ψₙ = normalize(L * ψ)
            # record jump time
            push!(jump_t_list, t)
            # record jump_operator
            push!(jump_L_list, L |> tag)
            # record time step for the ansatz
            push!(t_list, t + dt)
            # set new reference state
            set_ref!(A, ψₙ)
            # reset Ansatz
            reset_θ!(A)
            # record the new ansatz parameters after jump
            push!(θ_list, A.θ |> copy)
            # reset jump recorder
            Γ = 0
            q = rand(Float64)
        end
        save_state ? push!(u_list, A |> get_state) : nothing
        t += dt
        t_idx += 1
    end
    # record final ansatz operators and θ
    if (isempty(jump_t_list) || !isapprox(t, jump_t_list[end])) && !save_everystep
        push!(θ_list, A.θ |> copy)
        push!(A_list, [i |> tag for i in A.A])
    end
    # reset ansatz to the inital condition
    set_ref!(A, ref_init)
    reset!(A)
    # construct solution object
    AVQDSol(t_list, u_list, θ_list, A_list, jump_L_list, jump_t_list)
end

struct AVQDSol
    "time grid for the solution"
    t::Vector
    "state vector solution at each time grid point"
    u::Vector
    "variational parameters"
    θ::Vector
    "ansatz operators"
    A::Vector
    "jump operators"
    jump_L::Vector
    "jump time points"
    jump_t::Vector
    "exit status for the algorithm"
    status::Integer
    "norm of the state vector"
    norm::Vector
end

AVQDSol(t::Vector, u::Vector, θ::Vector, A::Vector, jump_L::Vector,
    jump_t::Vector; status::Integer=0, norm::Vector=[]) = AVQDSol(t, u, θ, A, jump_L, jump_t, status, norm)
