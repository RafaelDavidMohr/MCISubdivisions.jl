# --- mixed cells --- #

function find_nonzero_indices(V::Matrix{C}, inds::Vector{Int}, test_inds::Vector{Int},
                              V_rank::Int) where C

    res = Int[]
    Proj = project_along_linear_space(V[:, inds], V[:, test_inds], V_rank)
    for i in 1:size(Proj, 2)
        !iszero(Proj[:, i]) && push!(res, test_inds[i])
    end
    return res
end

function MixedCell(inds::Vector{Vector{Int}}, M::MCI)
    loc_inds = Vector{Int}[]
    for i in 1:length(inds)
        rem_inds = setdiff(indices(M), vcat(inds[1:i-1]...))
        rk = sum((length).(inds[1:i-1])) - (i - 1)
        loc_inds_i = isone(i) ? rem_inds : find_nonzero_indices(M.V, vcat(inds[1:i-1]...), rem_inds, rk)
        sort!(loc_inds_i)
        push!(loc_inds, loc_inds_i)
    end
    return MixedCell(inds, loc_inds)
end

function vol(m::MixedCell, A::Matrix{Int})
    res = 1
    for S in m.inds
        res *= lattice_volume(convex_hull(transpose(A[:, S])))
    end
    return res
end

function cayley_indices(m::MixedCell, j::Int)
    if j == length(m)
        return m.loc_inds[j]
    end
    res = setdiff(m.loc_inds[j], m.loc_inds[j+1])
    return sort(vcat(res, [first(m.inds[j+1])]))
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

    isempty(V) && return W

    F = parent(first(V))
    RV = Matrix(echelon_form(matrix(F, transpose(V))))

    pivot_cols = falses(size(RV, 2))
    pivot_rows = zeros(Int, size(RV, 2))
    @inbounds for i in 1:size(RV, 1)
        for j in 1:size(RV, 2)
            if !iszero(RV[i, j])
                pivot_cols[j] = true
                pivot_rows[j] = i
                break
            end
        end
    end

    non_pivot_count = size(RV, 2) - count(pivot_cols)
    non_pivot_inds = Vector{Int}(undef, non_pivot_count)
    idx = 1
    @inbounds for j in 1:size(RV, 2)
        if !pivot_cols[j]
            non_pivot_inds[idx] = j
            idx += 1
        end
    end

    result = Matrix{C}(undef, size(V, 1) - V_rank, size(W, 2))
    red = Vector{C}(undef, size(W, 1))

    @inbounds for i in 1:size(W, 2)
        for j in 1:size(W, 1)
            red[j] = W[j, i]
        end

        for j in 1:size(W, 1)
            if !iszero(red[j]) && pivot_cols[j]
                pind = pivot_rows[j]
                for k in j+1:size(W,1)
                    red[k] = red[k] - red[j] * RV[pind, k]
                end
                red[j] = F(0)
            end
        end
        result[:, i] = red[non_pivot_inds]
    end

    return result
end

function linear_span(A::Matrix{C}, inds::Vector{Int}) where C
    a0 = A[:, first(inds)]
    L = Matrix{C}(undef, size(A, 1), length(inds) - 1)
    for (j, i) in enumerate(inds[2:end])
        L[:, j] = A[:, i] - a0
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

function LinearAlgebra.dot(v::Vector{Float64}, c::Hyperplane)
    return dot(v, c.cfs_fl)
end

# --- homotopy paths --- #
    
function first_intersection(p0::Vector{Int}, p1::Vector{Int},
                            hyperplanes::AbstractSet{Hyperplane},
                            last_h::Union{Nothing, Hyperplane})

    best_h = nothing

    for c in hyperplanes
        !(does_cross(p0, p1, c)) && continue
        if isnothing(last_h) || !lt_refined(p0, p1, last_h, c)
            if isnothing(best_h) || lt_refined(p0, p1, best_h, c)
                best_h = c
            end
        end
    end

    return best_h
end
    
# c1 < c2 in lexicographic order refined by p0
function lt_refined(p0::Vector{Int}, p1::Vector{Int}, c1::Hyperplane, c2::Hyperplane)

    m1_fl, m1_P = dot(p1, c1) .^ (-1)
    m2_fl, m2_P = dot(p1, c2) .^ (-1)

    d1_fl, d1_P = (m1_fl, m1_P) .* dot(p0, c1)
    d2_fl, d2_P = (m2_fl, m2_P) .* dot(p0, c2)
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

function does_cross(p0::Vector{Int}, p1::Vector{Int}, c::Hyperplane)
    p0c, _ = dot(p0, c)
    p1c, p1c_P = dot(p1, c)
    iszero(p1c_P) && return false
    return signbit(p0c) ⊻ signbit(p1c)
end

function crossing_val(p0::Vector{Int}, p1::Vector{Int}, c::Hyperplane)
    p0c, _ = dot(p0, c)
    p1c, _ = dot(p1, c)
    return p0c / (p0c - p1c)
end

# --- other helpers --- #

# TODO probably: sorted merge
function restrict(inds::Vector{Int}, restr::Vector{Int})
    res = Int[]
    for i in inds
        !(i in restr) && continue
        push!(res, i)
    end
    return res
end

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
