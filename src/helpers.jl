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

function nz_inds(c::Hyperplane)
    return findall(!iszero, c.cfs_P)
end

function partial_sum(inds::Vector{Int}, c::Hyperplane)
    return sum(c.cfs_fl[inds]), sum(c.cfs_P[inds])
end

function partial_sum_sign(inds::Vector{Int}, c::Hyperplane, sgn::Bool)
    res, resP = partial_sum(inds, c)
    return resP != 0 && (sgn ? res > 0 : res < 0) # TODO: check if 0 is allowed here
end

function LinearAlgebra.dot(v::Vector{Int}, c::Hyperplane)
    return dot(v, c.cfs_fl), dot(v, c.cfs_P)
end

# --- homotopy paths --- #
    
function first_intersection(p0::Vector{Int}, p1::Vector{Int},
                            hyperplanes::AbstractSet{Hyperplane},
                            last_h::Union{Nothing, Hyperplane})

    mlh_fl, mlh_P = isnothing(last_h) ? (Inf, nothing) : dot(p1, last_h)

    best_h = nothing
    mbh_fl, mbh_P = Inf, nothing

    for c in hyperplanes
        mc_fl, mc_P = dot(p1, c)
        iszero(mc_P) && continue
        if isnothing(last_h) || !lt_refined(p0, c, last_h, mlh_fl, mlh_P, mc_fl, mc_P)
            if isnothing(best_h) || lt_refined(p0, c, best_h, mc_fl, mc_P, mbh_fl, mbh_P)
                best_h = c
                mbh_fl = mc_fl
                mbh_P = mc_P
            end
        end
    end

    return best_h
end
    
# c1 < c2 in lexicographic order refined by d flipped
function lt_refined(d::Vector{Int}, c1::Hyperplane, c2::Hyperplane,
                    m1_fl::Float64, m1_P::FqFieldElem,
                    m2_fl::Float64, m2_P::FqFieldElem)

    d1_fl, d1_P = (m1_fl, m1_P) .* dot(d, c1)
    d2_fl, d2_P = (m2_fl, m2_P) .* dot(d, c2)
    if d1_P != d2_P
        return d1_fl > d2_fl
    end

    for i in 1:length(c1)
        c1e_P = m1_P * c1.cfs_P[i]
        c2e_P = m2_P * c2.cfs_P[i]
        if c1e_P != c2e_P
            return m1_fl * c1.cfs_fl[i] > m2_fl * c2.cfs_fl[i]
        end
    end

    return true # error check here?
end

# --- other helpers --- #

function forgetful_lift(A_size::Int, forget_inds::Vector{Int})
    d = Vector{Int}(undef, A_size)
    for i in 1:A_size
        if i in forget_inds
            d[i] = 0
        else
            d[i] = 1
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
