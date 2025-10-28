function getindex(c::SparseVec, inds...)
    rem = findall(i -> i in inds, c.inds)
    return SparseVec(cfs[rem], c.inds[rem])
end

function subvec(c::SparseVec, inds...)
    rem = findall(i -> i in inds, c.inds)
    return view(c.cfs, rem)
end

function densify(c::SparseVec{C}, dim::Int) where C

    v = zeros(C, dim)
    j = 1
    for i in 1:n
        if i == c.inds[j]
            v[i] = c.cfs[j]
            j += 1
        end
    end
    return v
end

function add_to_dict!(d::Dict{T, Vector{S}}, k::T, v::S)
    if haskey(d, k)
        push!(d[k], v)
    else
        d[k] = [v]
    end
end

function is_affinely_independent(A_ext::Matrix{Int64})
    return rank(matrix(QQ, A_ext)) == size(A, 2)
end

function rank(M::MCI, I)
    return rank(matrix(QQ, M.V[:, I]))
end
