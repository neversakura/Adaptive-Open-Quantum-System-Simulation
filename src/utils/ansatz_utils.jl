"""
count_letters(str::AbstractString)

Returns the number of letters in a given string `str`.

# Arguments
- `str::AbstractString`: The input string to count letters in.

# Returns
An integer representing the number of letters in the input string.

# Examples
```julia
julia> str = "Hello, world!";
julia> count_letters(str)
10
"""
function count_letters(str::AbstractString)
    letter_count = 0
    for c in str
        if isletter(c)
            letter_count += 1
        end
    end
    return letter_count
end

"""
count_params(Alist, t, jump_t)

Given a list of ansatz strings `Alist`, a vector of time values `t`, and a
vector of jump times `jump_t`, computes the number of parameters and the number
of controlled-not gates in each ansatz string `A` at the jump times.

# Arguments
- `Alist`: a list of ansatz strings of length `N`, where `N` is the number
  of time steps.
- `t`: a vector of length `N` representing the time values.
- `jump_t`: a vector of length `M` representing the jump times.

# Returns
- `np_list`: a list of length `N` containing the number of parameters in each
  ansatz string `A` at the jump times.
- `cn_list`: a list of length `N` containing the number of controlled-not
  gates in each ansatz string `A` at the jump times.
"""
function count_params(Alist, t, jump_t)
    jpos = [find_positions(t, x) for x in jump_t]
    np_list = []
    cn_list = []
    for (A, posl) in zip(Alist, jpos)
        np = 0
        cn = 0
        for pos in posl
            Alen = length(A[pos])
            if Alen > 0
                np += Alen
                cn += sum((x) -> x == 2, [count_letters(x) for x in A[pos]])
            end
        end
        if isempty(posl) || posl[end] != length(A) - 1
            np += length(A[end])
            cn += sum((x) -> x == 2, [count_letters(x) for x in A[end]])
        end
        push!(np_list, np)
        push!(cn_list, 2 * cn)
    end
    np_list, cn_list
end

"""
find_positions(A::Vector, B::Vector)

Returns an array representing the positions in `A` of each element in `B`.

# Arguments
- `A::Vector`: The input array to search for elements from `B`.
- `B::Vector`: The input array of elements to search for in `A`.

# Returns
An array which contains the positions of an element in `B` in `A`.
"""
function find_positions(A::Vector, B::Vector)
    positions = []
    for b in B
        ps = findall(x -> isapprox(x, b), A)
        if length(ps) > 1
            error("Duplicated jump time.")
        end
        push!(positions, ps[1])
    end
    return positions
end

# =========== ansatz writting utility ==============
# writting parameters into CSV files
function write_arr_str(params::Vector{<:Number})
    if isempty(params)
        return "0"
    end
    res = reduce(*, [string(p) * " " for p in params])
    res[1:end-1]
end

function write_arr_str(params::Vector)
    if isempty(params)
        return "none"
    end
    res = reduce(*, [p * " " for p in params])
    res[1:end-1]
end

function write_jump_op_list(t_list, jump_t, jump_L)
    res = []
    for t in t_list
        idx = searchsortedlast(jump_t, t)

        if idx < 1
            push!(res, "none")
        else
            if jump_t[idx] ≈ t
                push!(res, jump_L[idx])
            else
                push!(res, "none")
            end
        end
    end
    res
end