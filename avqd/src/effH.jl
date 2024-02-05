struct EffectiveHamiltonian
    "Hermitian part"
    Hh::OpenQuantumBase.AbstractHamiltonian
    "antiHermitinan part"
    Ha::Matrix
    "List of L"
    L::Vector
    "List of L†L"
    LdL::Vector
end

function EffectiveHamiltonian(funcs, mats::Vector{<:AbstractMatrix}, γlist::Vector, Llist::Vector{<:TagOperator})
    H = Hamiltonian(funcs, mats, unit=:ħ, dimensionless_time=false)
    LdL = [g*L'*L for (g, L) in zip(γlist, Llist)]
    EffectiveHamiltonian(H, 0.5*sum(LdL), Llist, LdL)
end

function EffectiveHamiltonian(funcs, mats::Vector{<:AbstractMatrix}, γ::Number, Llist::Vector{<:TagOperator})
    H = Hamiltonian(funcs, mats, unit=:ħ, dimensionless_time=false)
    LdL = [γ*L'*L for L in Llist]
    EffectiveHamiltonian(H, 0.5*sum(LdL), Llist, LdL)
end

function EffectiveHamiltonian(funcs, mats::Vector{<:AbstractMatrix}, γlist::Vector{<:Number}, Llist::Vector{Vector{T}}) where T<:TagOperator
    H = Hamiltonian(funcs, mats, unit=:ħ, dimensionless_time=false)
    LdL = vcat([[γ*L'*L for L in LL] for (γ, LL) in zip(γlist, Llist)]...)
    EffectiveHamiltonian(H, 0.5*sum(LdL), vcat(Llist...), LdL)
end

herm(H::EffectiveHamiltonian, t) = H.Hh(t)
antiherm(H::EffectiveHamiltonian) = H.Ha
get_LdL(H::EffectiveHamiltonian) = H.LdL
get_L(H::EffectiveHamiltonian) = H.L

function vectorize_comm(A)
    iden = one(A)
    iden⊗A - transpose(A)⊗iden
end

struct VectorizedEffectiveHamiltonian
    "Hermitian part"
    Hh
    "antiHermitinan part"
    Ha::Matrix
end

function VectorizedEffectiveHamiltonian(funcs, mats::Vector{<:AbstractMatrix}, γ::Number, Llist::Vector{<:TagOperator})
    iden = one(mats[1])
    d = size(mats[1], 1)
    vec_mats = [vectorize_comm(m) for m in mats]
    res = zeros(ComplexF64, d^2, d^2)
    for L in Llist
        Ł = L'*L
        res = res - γ*(conj(L)⊗L-(iden⊗Ł+transpose(Ł)⊗iden)/2)
    end
    H = Hamiltonian(funcs, vec_mats, unit=:ħ, dimensionless_time=false)
    VectorizedEffectiveHamiltonian(H, res)
end

function VectorizedEffectiveHamiltonian(funcs, mats::Vector{<:AbstractMatrix}, γlist::Vector{<:Number}, Llist::Vector{Vector{T}}) where T<:TagOperator
    iden = one(mats[1])
    d = size(mats[1], 1)
    vec_mats = [vectorize_comm(m) for m in mats]
    res = zeros(ComplexF64, d^2, d^2)
    for (γ, LL) in zip(γlist, Llist)
        for L in LL
            Ł = L'*L
            res = res - γ*(conj(L)⊗L-(iden⊗Ł+transpose(Ł)⊗iden)/2)
        end
    end
    H = Hamiltonian(funcs, vec_mats, unit=:ħ, dimensionless_time=false)
    VectorizedEffectiveHamiltonian(H, res)
end

herm(H::VectorizedEffectiveHamiltonian, t) = H.Hh(t)
antiherm(H::VectorizedEffectiveHamiltonian) = H.Ha
