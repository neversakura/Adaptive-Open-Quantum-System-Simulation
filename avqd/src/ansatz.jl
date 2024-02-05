mutable struct Ansatz{relative_rcut}
    "Ansatz angles"
    θ
    "Current ansatz"
    A
    "Current state"
    state
    "Reference state"
    ref
    "Operator pools"
    pool
    "Maximum allowed McLachlan distance"
    rcut
    "Minimum allowed McLachlan distance improvement"
    relrcut
    "Maximum allowed ∂|ψ⟩/∂θᵢ"
    pthcut
    "Number of qubits"
    nqbit
end

function Ansatz(u0::Vector; θ=Float64[], ansatz=[], rcut=0.001, pthcut=10, pool="all2", vectorize=false, relrcut=nothing)
    nqbit = length(u0) |> log2 |> Int
    pool_nqbit = nqbit
    if vectorize
        tmp = u0 * u0'
        u0 = tmp[:]
        pool_nqbit = 2 * nqbit
    end
    pool = build_pool(pool_nqbit, type=pool)

    relative_rcut = relrcut === nothing ? false : true
    Ansatz{relative_rcut}(θ, ansatz, u0, copy(u0), pool, rcut, relrcut, pthcut, nqbit)
end

get_state(A::Ansatz) = A.state
pool_len(A::Ansatz) = length(A.pool)
pool_tag(A::Ansatz) = [tag(i) for i in A.pool]
get_ref(A::Ansatz) = A.ref
get_newest_A(A::Ansatz) = isempty(A.A) ? nothing : A.A[end]
update_θ!(A::Ansatz, dθ, dt) = A.θ = A.θ + dθ * dt
reset_θ!(A::Ansatz) = fill!(A.θ, 0)
set_ref!(A::Ansatz, ref) = A.ref = ref

function set_pool_to_ansatz!(A::Ansatz)
    A.A = copy(A.pool)
    A.θ = zeros(length(A.pool))
end

function add_A!(A::Ansatz, P::AnsatzOperatorBase)
    push!(A.A, P)
    push!(A.θ, 0)
end

function reset!(A::Ansatz)
    A.θ = Float64[]
    A.A = []
    update_state!(A)
end

function update_state!(A::Ansatz)
    if !isempty(A.A)
        state = A.ref |> copy
        for i in 1:length(A.A)
            state = aexp(A.A[i], A.θ[i]) * state
        end
        A.state = state
    else
        A.state = A.ref
    end
end

function partial_theta(A::Ansatz)
    len = length(A.A)
    ψ = copy(A.ref)
    res = []
    for i = 1:len
        ψ = aexp(A.A[i], A.θ[i]) * ψ
        ψt = -0.5im * lmul(A.A[i], ψ)
        for j = i+1:len
            ψt = aexp(A.A[j], A.θ[j]) * ψt
        end
        push!(res, ψt)
    end
    res
end

function one_step!(A::Ansatz{true}, He, Ha, dt)
    relrcut = A.relrcut
    ψ = A |> get_state
    dψ = partial_theta(A)
    M = build_m(ψ, dψ)
    V = build_v(ψ, dψ, He, Ha)

    dθ, vmv = lin_solve(M, V)

    vmvMax = vmv
    Mtmp = M
    Vtmp = V
    dθtmp = dθ
    opTmp = nothing
    dψₐTmp = nothing

    add_flag = true

    while add_flag
        tagTmp = nothing
        for op in A.pool
            if (get_newest_A(A) |> tag) == tag(op)
                continue
            end

            dψₐ = -0.5im * lmul(op, ψ)
            Mop = update_m(M, dψₐ, ψ, dψ)
            Vop = update_v(V, dψₐ, ψ, He, Ha)
            dθop, vmvOp = lin_solve(Mop, Vop)

            if vmvOp > vmvMax
                # adding op decrease the distance
                Mtmp = Mop
                Vtmp = Vop
                vmvMax = vmvOp
                dθtmp = dθop
                opTmp = op
                tagTmp = tag(op)
                dψₐTmp = dψₐ
            end
        end
        add_flag = vmvMax - vmv < relrcut ? false : true
        if tagTmp !== nothing && add_flag
            #@info "Add operator to ansatz" tagTmp
            add_A!(A, opTmp)
            vmv = vmvMax
            M = Mtmp
            V = Vtmp
            dθ = dθtmp
            push!(dψ, dψₐTmp)
        end
    end
    update_θ!(A, dθ, dt)
    update_state!(A)
end

function one_step!(A::Ansatz{false}, He, Ha, dt)
    ψ = A |> get_state
    dψ = partial_theta(A)
    M = build_m(ψ, dψ)
    V = build_v(ψ, dψ, He, Ha)

    dθ, vmv = lin_solve(M, V)
    varHe = 2 * real(ψ' * He * He * ψ - (ψ' * He * ψ)^2)

    vmvMax = vmv
    Mtmp = M
    Vtmp = V
    dθtmp = dθ
    opTmp = nothing
    tagTmp = nothing
    dψₐTmp = nothing
    while varHe - vmv > A.rcut
        tagTmp = nothing
        for op in A.pool
            if (get_newest_A(A) |> tag) == tag(op)
                continue
            end

            dψₐ = -0.5im * lmul(op, ψ)
            Mop = update_m(M, dψₐ, ψ, dψ)
            Vop = update_v(V, dψₐ, ψ, He, Ha)
            dθop, vmvOp = lin_solve(Mop, Vop)
            if vmvOp > vmvMax
                # adding op decrease the distance
                Mtmp = Mop
                Vtmp = Vop
                dθtmp = dθop
                vmvMax = vmvOp
                opTmp = op
                tagTmp = tag(op)
                dψₐTmp = dψₐ
            end
        end

        tagTmp === nothing ? throw("No operator in ansatz pool can improve the McLachlan distance below threshold.") : nothing


        @debug "Add operator to ansatz" tagTmp
        add_A!(A, opTmp)
        vmv = vmvMax
        M = Mtmp
        V = Vtmp
        dθ = dθtmp
        push!(dψ, dψₐTmp)
        @debug "McLachlan distance decreases to" varHe - vmv + bound
    end
    update_θ!(A, dθ, dt)
    update_state!(A)
end

function one_step_non_adp!(A::Ansatz, He, Ha, dt)
    ψ = A |> get_state
    dψ = partial_theta(A)
    M = build_m(ψ, dψ)
    V = build_v(ψ, dψ, He, Ha)

    dθ, = lin_solve(M, V)

    update_θ!(A, dθ, dt)
    update_state!(A)
end

function build_m(ψ, dψ)
    l = length(dψ)
    res = zeros(l, l)
    for μ in 1:l
        for ν in 1:l
            res[μ, ν] = 2 * real(dψ[μ]' * dψ[ν] + ψ' * dψ[μ] * ψ' * dψ[ν])
        end
    end
    res
end

function update_m(m, dψₐ, ψ, dψ)
    l = size(m, 1)
    mp = zeros(l + 1, l + 1)
    mp[l+1, l+1] = 2 * real(dψₐ' * dψₐ + ψ' * dψₐ * ψ' * dψₐ)
    mp[1:l, 1:l] = m
    for μ in 1:l
        temp = 2 * real(dψ[μ]' * dψₐ + ψ' * dψ[μ] * ψ' * dψₐ)
        mp[μ, l+1] = temp
        mp[l+1, μ] = temp
    end
    mp
end

function build_v(ψ, dψ, He, Ha)
    L = -1im * (He - 1im * Ha)
    l = length(dψ)
    res = zeros(l)
    for μ = 1:l
        res[μ] = 2 * real(dψ[μ]' * L * ψ + dψ[μ]' * ψ * ψ' * L' * ψ)
    end
    res
end

function update_v(v, dψₐ, ψ, He, Ha)
    L = -1im * (He - 1im * Ha)
    l = length(v)
    vp = zeros(l + 1)
    vp[1:l] = v
    vp[l+1] = 2 * real(dψₐ' * L * ψ + dψₐ' * ψ * ψ' * L' * ψ)
    vp
end

"""
    build_pool(nqbit::Integer; type::String="all2")

Build a `nqbit` operator pool specified by `type`. If `nqbit` is greater than 4, the function will automatically switch to sparse matrix.

# Arguments
- `nqbit::Integer`: total number of qubit
- `type::String="all2"`: type of ansatz pool to construct. There are three options: "all2": a pool with every single-qubit and two-qubit Pauli operator; "neighbor": a pool with every single-qubit operator and all the two-qubit Pauli operator on the nearest neighboring qubits of a 1-D chain; "ASC": a pool with every single-qubit operator and all the two-qubit ZZ operator on the nearest neighboring qubits of a 1-D chain.
"""
function build_pool(nqbit::Integer; type::String="all2")
    ltype = lowercase(type)
    if ltype == "all2"
        build_pool_all(nqbit)
    elseif ltype == "neighbor"
        build_pool_neighbor(nqbit)
    elseif ltype == "asc"
        build_pool_asc(nqbit)
    elseif ltype == "vectorized_neighbor"
        build_pool_vectorized_neighbor(nqbit)
    else
        throw(ArgumentError("Pool type $type not supported."))
    end
end

"""
    build_pool_all(nqbit::Int)

Build an operator pool with every single-qubit Pauli operator and every two-qubit Pauli operator.
"""
function build_pool_all(nqbit::Int)
    pauliStr = ["X", "Z", "Y"]
    res = []
    for order = 1:2
        for idx in combinations(1:nqbit, order)
            for op in Iterators.product(repeat([pauliStr], order)...)
                push!(res, PauliOperator([i for i in op], idx, 1, nqbit))
            end
        end
    end
    res
end

"""
    build_pool_neighbor(nqbit::Int)

Build an operator pool with every single qubit Pauli operator and all the two-qubit Pauli operator on the nearest neighboring qubits of a 1-D chain.
"""
function build_pool_neighbor(nqbit::Int)
    pauliStr = ["X", "Z", "Y"]
    res = []
    for idx = 1:nqbit
        for op in pauliStr
            push!(res, PauliOperator([op], [idx], 1, nqbit))
        end
        if idx + 1 <= nqbit
            for op in Iterators.product(repeat([pauliStr], 2)...)
                push!(res, PauliOperator([i for i in op], [idx, idx + 1], 1, nqbit))
            end
        end
    end
    res
end

"""
    build_pool_asc(nqbit::Int)

Build an operator pool with every single-qubit Pauli operator and every two-qubit ZZ operator on the nearest neighboring qubits of a 1-D chain.
"""
function build_pool_asc(nqbit::Int)
    pauliStr = ["X", "Z", "Y"]
    res = []
    for idx = 1:nqbit
        for op in pauliStr
            push!(res, PauliOperator([op], [idx], 1, nqbit))
        end
        if idx + 1 <= nqbit
            push!(res, PauliOperator(["Z", "Z"], [idx, idx + 1], 1, nqbit))
        end
    end
    res
end

function build_pool_vectorized_neighbor(nqbit::Int)
    nvqbit = 2 * nqbit
    pauliStr = ["X", "Z", "Y"]
    res = []
    for idx = 1:nvqbit
        for op in pauliStr
            push!(res, PauliOperator([op], [idx], 1, nqbit))
        end
        if (idx + 1 <= nvqbit) && (idx != nqbit)
            for op in Iterators.product(repeat([pauliStr], 2)...)
                push!(res, PauliOperator([i for i in op], [idx, idx + 1], 1, nqbit))
            end
        end
        if idx <= nqbit
            for op in Iterators.product(repeat([pauliStr], 2)...)
                push!(res, PauliOperator([i for i in op], [idx, idx + nqbit], 1, nqbit))
            end
        end
    end
    res
end

function lin_solve(m::AbstractMatrix{Float64}, v::AbstractVector{Float64}, λ::Float64=1e-4)
    A = m'm
    xλ = try
        cholesky!(Hermitian(zot(A, λ^2.0))) \ (m' * v)
    catch
        zot(A, λ^2.0) \ (m' * v)
    end
    xλ, 2 * v' * xλ - xλ' * m * xλ
end

function zot(A::AbstractMatrix{<:AbstractFloat}, λ::AbstractFloat)
    a = deepcopy(A)
    n = size(A'A, 1)
    for i = 1:n
        @inbounds a[i, i] += λ
    end
    return a
end
