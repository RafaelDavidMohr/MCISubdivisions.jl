# --- matrix functions --- #

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
    A_aff = vcat(A[:, inds], transpose([F(1) for _ in 1:length(inds)]))
    return rank(matrix(F, A_aff)) == length(inds)
end

# --- helper for our Circuit data structure --- #

function partial_sum(inds::Vector{Int}, c::Circuit)
    res = 0.0
    for (i, ind) in enumerate(c.inds)
        if ind in inds
            res += c.cfs_fl[i]
        end
    end
    return res
end

function LinearAlgebra.dot(v::Vector{Float64}, c::Circuit)
    res = 0.0
    i = 1
    for j in 1:length(v)
        i > length(c.inds) && break
        if j == c.inds[i]
            res += v[j]*c.cfs_fl[i]
            i += 1
        end
    end
    return res
end

# first_intersection(p0, p1, hyperplanes)

# given two points p0, p1 and hyperplanes hyperplanes = [ ([a1, a2, ... , an]), ... , )
# return the hyperplane that p0 + t(p1 - p0), t = [0,1], intersects the first, the value of t on this intersection 
# and the point of the intersection 
    
function first_intersection(p0::Vector{Float64},
                            p1::Vector{Float64},
                            hyperplanes::Vector{Circuit}) 

    d = p1 - p0

    best_t = Inf
    best_h = nothing

    for c in hyperplanes
        denom = dot(d, c)
        @assert abs(denom) < 1e-12 "Path not generic enough!"
        t = -(dot(p0, c)) / denom
        if 0 < t <= 1 && t < best_t # t = 0 explicitly excluded
            best_t = t
            best_h = c
        end
    end

    return best_t, best_h
end

# to compute the total length of the given piecewise linear path
function path_length(path::HomotopyPath)
    points = path.points
    total = 0.0
    for i in 1:(length(points)-1)
        total += norm(points[i+1] - points[i])
    end
    return total
end

# to compute the piecewise linear path on n = size(p0) points
function piecewice_linear_path(p0::Vector{QQFieldElem},
                               p1::Vector{QQFieldElem}) 

    n = length(p0)
    eps = 1e-8
    points = Vector{Float64}[]

    push!(points, p0)

    prev = p0
    for _ in 1:n
        t = rand()  # random number in [0,1]
        shifts = (2 .* rand(n) .- 1) .* eps #random numbers in [-eps,eps]
        # Next point lies between prev and p1
        next_point = prev + t * (p1 - prev) + shifts 
        push!(points, next_point)
        prev = next_point
    end

    push!(points, p1)

    return HomotopyPath(points)
end

# to compute the point of the first intersection of piecewise linear path with given hyperplanes
function first_intersection_with_path!(path::HomotopyPath,
                                       hyperplanes::AbstractSet{Circuit}) 

    n = length(path)

    t_int, h_int = Inf, nothing
    while isnothing(h_int) && !is_completed(path)
        p1, p2 = path[1], path[2]
        t_int, h_int = first_intersection(p1, p2, hyperplanes)
        isnothing(h_int) && popfirst!(path.points)
    end

    return t_int, h_int
end

# --- small helpers --- #

function add_to_dict!(d::Dict{T, Set{S}}, k::T, v::S) where {T, S}
    if haskey(d, k)
        push!(d[k], v)
    else
        d[k] = Set([v])
    end
end

# random point on line segment between p1 and p2 
# with upper bound on parameter t
function random_path_point(p1::Vector{Float64},
                           p2::Vector{Float64},
                           upb_t::Float64)

    denom = ceil(Int, 1000/upb_t)
    a = rand(1:denom-1)
    t = a/denom
    p_fl = (1 - t)*p1 + t*p2
    return (x -> QQ(rationalize(Int32, p_fl))).(p_fl)
end
