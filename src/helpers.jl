function add_to_dict!(d::Dict{T, Set{S}}, k::T, v::S) where {T, S}
    if haskey(d, k)
        push!(d[k], v)
    else
        d[k] = Set([v])
    end
end

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


# first_intersection(p0, p1, hyperplanes)

# given two points p0, p1 and hyperplanes hyperplanes = [ ([a1, a2, ... , an]), ... , )
# return the hyperplane that p0 + t(p1 - p0), t = [0,1], intersects the first, the value of t on this intersection 
# and the point of the intersection 
    
function first_intersection(p0::AbstractVector, p1::AbstractVector, hyperplanes::Vector{<:AbstractVector}) 

    d = float(p1) .- float(p0)
    fp0 = float(p0)

    best_t = Inf
    min_t = Inf
    best_i = nothing

    for (i, a) in enumerate(hyperplanes)

        a = float(a)
        denom = dot(a, d)
        if abs(denom) < 1e-12
            error("Path is not generic enought!")
        end
        # t = (-a p0) / (a (p1 - p0))
        t = -(dot(a, fp0)) / denom
        if 0 <= t <= 1 && t < best_t
            best_t = t
            best_i = i
        end
    end

    if isnothing(best_i)
        @info "no intersection with given hyperplanes"
        return nothing
    else
        min_t = -(dot(hyperplanes[best_i], p0)) / dot(hyperplanes[best_i], p1 .- p0)
        x = p0 + min_t * d
        return (best_i, min_t, x)
    end
end

# to compute the total length of the given piecewise linear path
function path_length(points::Vector{<:AbstractVector})
    total = 0.0
    for i in 1:(length(points)-1)
        total += norm(points[i+1] .- points[i])
    end
    return total
end

# to compute the piecewise linear path on n = size(p0) points
function piecewice_linear_path(p0::AbstractVector, p1::AbstractVector) 

    n = length(p0)
    eps = 1e-8
    points = Vector{AbstractVector}()

    push!(points, p0)

    prev = p0
    for _ in 1:n
        t = rand()  # random number in [0,1]
        shifts = (2 .* rand(n) .- 1) .* eps #random numbers in [-eps,eps]
        # Next point lies between prev and p1
        next_point = prev .+ t .* (p1 .- prev) .+ shifts 
        push!(points, next_point)
        prev = next_point
    end

    push!(points, p1)

return points
end

# to compute the point of the first intersection of piecewise linear path with given hyperplanes
function first_intersection_linear_path(path::Vector{<:AbstractVector}, hyperplanes::Vector{<:AbstractVector}) 

    n = length(path)

    for i in 1:(n-1)
        p1, p2 = path[i], path[i+1]
        result = first_intersection(p1, p2, hyperplanes)
        if result !== nothing
            return result
        end
    end

    return nothing
end
