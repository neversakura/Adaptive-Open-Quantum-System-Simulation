abstract type AnsatzOperatorBase end

tag(::Nothing) = nothing
tag(A::AnsatzOperatorBase) = A.tag
nqubit(A::AnsatzOperatorBase) = A.nqbit
(A::AnsatzOperatorBase)() = A.mat

lmul(A::AnsatzOperatorBase, ψ) = A.mat * ψ
aexp(A::AnsatzOperatorBase, θ) = exp(-1im*θ*A.mat/2)

Base.size(A::AnsatzOperatorBase) = size(A.mat)
Base.size(A::AnsatzOperatorBase, d) = size(A.mat, d)
Base.:*(a::Number, p::AnsatzOperatorBase) = a * p.mat
Base.:*(p::AnsatzOperatorBase, a::Number) = p.mat * a
Base.:*(v::AbstractVecOrMat, p::AnsatzOperatorBase) = v * p.mat
Base.:*(p::AnsatzOperatorBase, v::AbstractVecOrMat) = p.mat * v
Base.:*(a::AnsatzOperatorBase, b::AnsatzOperatorBase) = a.mat * b.mat
Base.adjoint(a::AnsatzOperatorBase) = a.mat'
Base.conj(A::AnsatzOperatorBase) = conj(A.mat)
Base.kron(A::AnsatzOperatorBase, b::AbstractMatrix) = kron(A.mat, b)
Base.kron(b::AbstractMatrix, A::AnsatzOperatorBase) = kron(b, A.mat)

struct PauliOperator <: AnsatzOperatorBase
    mat
    tag
    nqbit
end

aexp(P::PauliOperator, θ) = cos(θ/2) * I - 1.0im * sin(θ/2) * P.mat

function PauliOperator(ops::Vector, idx::Vector, w::Number, nqbit::Int)
    sortedRep = sort([[o, i] for (o, i) in zip(ops, idx)], by=(x)->x[2])
    tag = ""
    for (o, i) in sortedRep
        tag *= o * string(i)
    end
    sp = nqbit > 4 ? true : false
    # temporary solution: use sparse matrix for more than 4 qubits
    mat = single_clause(ops, idx, w, nqbit, sp=sp)
    PauliOperator(mat, tag, nqbit)
end

struct AnsatzOperator <: AnsatzOperatorBase
    mat
    tag
    nqbit
end

function pauli_str_split(s)
    qsplit = r"([XYZ])([0-9]+)"
    res = []
    for qs in  eachmatch(qsplit, s)
        push!(res, [qs.captures[1], qs.captures[2]])
    end
    res
end

function vectorize(p_str, nqbit)
	p_str_list = pauli_str_split(p_str)
	first_term = []
	sign = -1
	for p in p_str_list
		pos = parse(Int, p[2]) + nqbit
		push!(first_term, [p[1], string(pos)])
		sign = p[1] == "Y" ? -1*sign : sign
	end
	sign_sym = sign > 0 ? "+" : "-"
	reduce(*, Iterators.flatten(first_term)) * sign_sym * p_str
end

function vectorize(P::PauliOperator)
    p_tag = P |> tag
    nqbit = P |> nqubit
    n_tag = vectorize(p_tag, nqbit)
    mat = vectorize_comm(P.mat)
    AnsatzOperator(mat, n_tag, 2^(2*nqbit))
end

function lind_vectorize(p_str, nqbit)
	p_str_list = pauli_str_split(p_str)
	second_term = []
	for p in p_str_list
		pos = parse(Int, p[2]) + nqbit
		push!(second_term, [p[1], string(pos)])
	end
	p_str * reduce(*, Iterators.flatten(second_term))
end

struct TagOperator <: AnsatzOperatorBase
    mat
    tag::String
    nqbit::Integer
end

function count_ops(ops::Vector)
    num_ops = [pauli_str_split(op)|>length for op in ops]
    res = zeros(Int, maximum(num_ops))
    for n in num_ops
        res[n] += 1
    end
    res
end
