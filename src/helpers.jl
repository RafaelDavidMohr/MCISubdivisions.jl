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

# --- matrix functions --- #

function reduce_mod_rand_prime(V::Matrix{QQFieldElem})
    p = Hecke.rand_bits_prime(ZZ, 31)
    F = GF(p)
    return F, [GF(numerator(x)) * GF(denominator(x))^(-1) for x in V]
end

function rank(V::Matrix{QQFieldElem}, I)
    F, Vp = reduce_mod_rand_prime(V[:, I])
    return rank(matrix(F, Vp))
end

# first_intersection(p0, p1, hyplanes)

# given two points p0, p1 and hyperplanes hyplanes = [ ([a1, a2, ... , an]), ... , )
# return the hyperplane that p0 + t(p1 - p0), t = [0,1], intersects the first, the value of t on this intersection 
# and the point of the intersection 
    
function first_intersection(p0, p1, hyplanes)

    d = p1 - p0

    best_t = Inf
    best_i = nothing

    for (i, a) in enumerate(hyplanes)
        denom = dot(a, d)
        if abs(denom) < 1e-12
            error("Path is not generic enought!")
        end
        # t = (-a p0) / (a (p1 - p0))
        t = -(dot(a, p0)) / denom
        if 0 <= t <= 1 && t < best_t
            best_t = t
            best_i = i
        end
    end

    if isnothing(best_i)
        @info "no intersection with given hyperplanes"
        return (nothing, nothing, nothing)
    else
        x = p0 + best_t * d
        return (best_i, best_t, x)
    end
end
