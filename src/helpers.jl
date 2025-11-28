# --- mixed cells --- #

function vol(m::MixedCell, A::Matrix{Int})
    res = 1
    for S in m.inds
        res *= lattice_volume(convex_hull(transpose(A[:, S])))
    end
    return res
end

# --- matrix --- #

function reduce_mod_rand_prime(V::Matrix{QQFieldElem})
    p = Hecke.rand_bits_prime(ZZ, 31)
    F = GF(p)
    return [F(numerator(x)) * F(denominator(x))^(-1) for x in V]
end

function reduce_mod_rand_prime(V::Matrix{Int})
    p = Hecke.rand_bits_prime(ZZ, 31)
    F = GF(p)
    return [F(x) for x in V]
end

function project_along_linear_space(V::Matrix{C}, W::Matrix{C}, V_rank::Int) where C

    VW = hcat(V, W)
    F = parent(first(V))
    R = echelon_form(matrix(F, VW), reduced = false)
    return Matrix(R[V_rank + 1:end, size(V, 2) + 1:end])
end

function linear_span(A::Matrix{C}, inds::Vector{Int}) where C
    a0 = A[:, first(inds)]
    L = Matrix{C}(undef, size(A, 1), 0)
    for i in inds[2:end]
        L = hcat(L, A[:, i] - a0)
    end
    return L
end

function is_affine_independent(A::Matrix{C}, inds::Vector{Int}) where C
    F = parent(first(A))
    L = linear_span(A, inds)
    return rank(matrix(F, L)) == length(inds) - 1
end

# --- circuits --- #

# potentially to optimize
function partial_sum(inds::Vector{Int}, c::Circuit)
    res = 0.0
    resP = parent(first(c.cfs_modP))(0)
    for (i, ind) in enumerate(c.inds)
        if ind in inds
            res += c.cfs_fl[i]
            resP += c.cfs_modP[i]
        end
    end
    return res, resP
end

function partial_sum_sign(inds::Vector{Int}, c::Circuit, sgn::Bool)
    res, resP = partial_sum(inds, c)
    return resP != 0 && (sgn ? res > 0 : res < 0)
end

function LinearAlgebra.dot(v::Vector{Float64}, c::Circuit)
    res = 0.0
    for (i, ind) in enumerate(c.inds)
        res += v[ind]*c.cfs_fl[i]
    end
    return res
end

# --- homotopy paths --- #
    
function first_intersection(p0::Vector{Float64},
                            p1::Vector{Float64},
                            t0::Float64,
                            hyperplanes::AbstractSet{Circuit},
                            h_excluded::Union{Nothing, Circuit}) 

    d = p1 - p0

    best_t = Inf
    best_h = nothing

    for c in hyperplanes
        c == h_excluded && continue
        denom = dot(d, c)
        @assert abs(denom) > 1e-12 "Path not generic enough!"
        t = -(dot(p0, c)) / denom
        if t0 < t <= 1 && t < best_t # t0 explicitly excluded
            best_t = t
            best_h = c
        end
    end

    return best_t, best_h
end

function path_length(path::HomotopyPath)
    points = path.points
    total = 0.0
    for i in 1:(length(points)-1)
        total += norm(points[i+1] - points[i])
    end
    return total
end

function piecewice_linear_path(p0::Vector{Float64},
                               p1::Vector{Float64},
                               nsegs=length(p0)::Int) 

    eps = 1e-2
    points = Vector{Float64}[]

    push!(points, p0)
    for i in 1:nsegs-1
        shifts = (2 .* rand(length(p0)) .- 1) .* eps # random numbers in [-eps,eps]
        pt = (p0 .+ (i/nsegs) .* (p1 - p0)) + shifts
        push!(points, pt)
    end
    push!(points, p1)

    return HomotopyPath(0.0, points)
end

function first_intersection_with_path!(path::HomotopyPath,
                                       hyperplanes::AbstractSet{Circuit},
                                       h_excluded::Union{Nothing, Circuit}) 

    n = length(path)

    h_int = nothing
    while isnothing(h_int) && !is_completed(path)
        p1, p2 = path[1], path[2]
        t_int, h_int = first_intersection(p1, p2, path.t_curr, hyperplanes, h_excluded)
        if isnothing(h_int)
            popfirst!(path.points)
            path.t_curr = 0.0
        else
            path.t_curr = t_int
        end
    end

    return h_int
end

# --- other helpers --- #

function forgetful_lift(A_size::Int, forget_inds::Vector{Int})
    eps = 1e-3
    d = Vector{Float64}(undef, A_size)
    w = rand(50:100, A_size)
    for i in 1:A_size
        if i in forget_inds
            d[i] = -100 + eps*w[i]
        else
            d[i] = eps*w[i]
        end
    end
    return d
end

function get_eci_data(F::Vector{<:MPolyRingElem})
    exps = unique(vcat([collect(exponents(f)) for f in F]...))
    A = hcat(exps...)
    CT = eltype(base_ring(parent(first(F))))
    V = [coeff(f, A[:, i]) for f in F, i in 1:size(A, 2)]
    return A, V
end

function rand_vec_ff(F::FqField, n::Int)
    return (F).(rand(0:characteristic(F)-1, n))
end

function delete_mixed_cell!(wd::WalkData, m::MixedCell)
    for c in keys(wd.walls)
        filter!(((m0, i, j),) -> m != m0, wd.walls[c])
        if isempty(wd.walls[c])
            delete!(wd.walls, c)
        end
    end
end

function add_to_dict!(d::Dict{T, Set{S}}, k::T, v::S) where {T, S}
    if haskey(d, k)
        push!(d[k], v)
    else
        d[k] = Set([v])
    end
end

function compose_as_maps(m1::Vector{Int}, m2::Vector{Int})
    result = similar(m2)
    for (j, i) in enumerate(m2)
        result[j] = m1[i]
    end
    return result
end
