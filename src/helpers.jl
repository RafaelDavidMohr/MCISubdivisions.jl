# --- circuit computation --- #

function exchange_circuits(Mloc::RelativeMCI, known_circuit_inds::Vector{Int})
    F = prime_field_V(Mloc)
    inds = setdiff(1:size(Mloc.V_rel, 2), known_circuit_inds)
    if rank(matrix(F, Mloc.V_rel[:, inds])) < length(inds)
        return [inds]
    end
    result = Set{Vector{Int}}()
    for i in known_circuit_inds
        mat = matrix(F, Mloc.V_rel[:, 1:end .!= i])
        K = kernel(mat, side = :right)
        @assert size(K, 2) == 1
        new_c = Int[]
        for j in 1:size(K, 1)
            if !iszero(K[j, 1])
                j >= i ? push!(new_c, j + 1) : push!(new_c, j)
            end
        end
        push!(result, new_c)
    end
    return [c for c in collect(result)]
end

# --- mixed cells --- #

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

function find_nonzero_indices(V::Matrix{C}, inds::Vector{Int}, test_inds::Vector{Int},
                              V_rank::Int) where C

    res = Int[]
    F = parent(first(V))
    for i in test_inds
        if rank(matrix(F, V[:, vcat(inds, [i])])) == V_rank + 1
            push!(res, i)
        end
    end
    return res
end

function vol(m::MixedCell, A::Matrix{Int})
    mat = hcat([linear_span(A, S) for S in m.inds]...)
    lu = LinearAlgebra.lu(mat)
    return round(Int, abs(LinearAlgebra.det(lu)))
end

function real_root_count(m::MixedCell, A::Matrix{Int}, V::Matrix{QQFieldElem})
    F = GF(2)
    n = size(A, 1)
    ml = sum((length).(m.inds))
    ev = matrix(F, vcat(A[:, vcat(m.inds...)], ones(Int, 1, ml)))
    k = permutedims(Matrix(kernel(ev)))
    enc = size(k, 2)
    sol = try
        solve(ev, sgn(m, V))
    catch ArgumentError
        return 0
    end
    size(k, 2) == 0 && return 1
    result_set = Set{Vector{FqFieldElem}}()
    for i in 0:(2^enc -1)
        new_sol = sol + k*digits(i, base=2, pad=enc)
        push!(result_set, new_sol[1:n])
    end
    return length(result_set)
end

function sgn(m::MixedCell, V::Matrix{QQFieldElem})
    F = GF(2)
    sgns = eltype(F)[]
    for (i, S) in enumerate(m.inds)
        for j in S
            V_submat = V[:, vcat([p == i ? [l for l in S if l != j] : m.inds[p][2:end] for p in 1:length(m.inds)]...)]
            dt = det(matrix(QQ, V_submat))
            dt > 0 ? push!(sgns, F(1)) : push!(sgns, F(0))
        end
    end
    return sgns
end

function real_root_count(ms::Vector{MixedCell}, A::Matrix{Int}, V::Matrix{QQFieldElem})
    return sum([real_root_count(m, A, V) for m in ms])
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
    return c.nzinds
end

function partial_sum(inds::Vector{Int}, c::Hyperplane)
    return sum(c.cfs_fl[inds]), sum(c.cfs_P[inds])
end

function partial_sum_sign(inds::Vector{Int}, c::Hyperplane, sgn::Bool)
    res, resP = partial_sum(inds, c)
    return resP != 0 && (sgn ? res > 0 : res < 0) # TODO: check if 0 is allowed here
end

function LinearAlgebra.dot(v::Vector{Int}, c::Hyperplane)
    F = parent(first(c.cfs_P))
    res_P = F(0)
    res_fl = 0.0
    for i in nz_inds(c)
        res_P += v[i] * c.cfs_P[i]
        res_fl += v[i] * c.cfs_fl[i]
    end
    return res_fl, res_P
end

function LinearAlgebra.dot(v::Vector{Float64}, c::Hyperplane)
    res_fl = 0.0
    for i in nz_inds(c)
        res_fl += v[i] * c.cfs_fl[i]
    end
    return res_fl
end

function LinearAlgebra.dot(l::DualVector, c::Hyperplane)
    return DualNumber(dot(l.r, c)..., dot(l.eps, c)...)
end

# --- homotopy paths --- #
    
function first_intersection!(l0::DualVector, l1::DualVector,
                             hyperplanes::AbstractSet{Hyperplane},
                             last_h::Union{Nothing, Hyperplane},
                             wd::WalkData)

    F = prime_field_A(wd.M)
    one_dual = dual_one(F)

    last_t = if isnothing(last_h)
        dual_zero(F)
    else
        crossing_val(l0, l1, last_h)
    end

    best_t = nothing
    best_h = nothing

    for c in hyperplanes
        c in wd.nocross && continue
        if c == last_h
            push!(wd.nocross, c)
            continue
        end

        t_c = crossing_val(l0, l1, c)
        # check if crossing between last t and 1
        if lt_dual(t_c, last_t) || lt_dual(one_dual, t_c)
            push!(wd.nocross, c)
            continue
        end

        if isnothing(best_h) || lt_dual(t_c, best_t)
            best_t = t_c
            best_h = c
        end
    end

    return best_h, best_t
end

function crossing_val(l0::DualVector, l1::DualVector, c::Hyperplane)
    l0d = dot(l0, c)
    l1d = dot(l1, c)
    denom = l0d - l1d

    F = prime_field(l0d)
    try
        if iszero(l0d.r_P) && iszero(l1d.r_P)
            return DualNumber(l0d.eps_fl * denom.eps_fl^(-1),
                              l0d.eps_P * denom.eps_P^(-1), 0.0, F(0))
        elseif iszero(denom.r_P)
            return DualNumber(2.0, F(2), 0.0, F(0))
        else
            return l0d * inv(denom)
        end
    catch e
        println(l0)
        println(l1)
        println(c.cfs_P)
        println(denom)
        println(characteristic(F))
        rethrow(e)
    end
end

# --- other helpers --- #

function to_mat_dual(a::Tuple{T, T}) where T
    return [a[1] a[2]; zero(a[1]) a[1]]
end

function to_dual_mat(a::Matrix{T}) where T
    return a[1, 1], a[1, 2]
end

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
    return DualVector(d)
end

function get_support(F::Vector{<:MPolyRingElem})
    exps = unique(vcat([collect(exponents(f)) for f in F]...))
    return hcat(exps...)
end

function get_eci_data(F::Vector{<:MPolyRingElem})
    A = get_support(F)
    CT = eltype(base_ring(parent(first(F))))
    V = [coeff(f, A[:, i]) for f in F, i in 1:size(A, 2)]
    return A, V
end

function rand_arr_ff(F::FqField, dims...)
    return (F).(rand(0:characteristic(F)-1, dims...))
end

function delete_mixed_cell!(wd::WalkData, m::MixedCell)
    for c in keys(wd.walls)
        filter!(((m0, i, j),) -> m != m0, wd.walls[c])
        if isempty(wd.walls[c])
            delete!(wd.walls, c)
        end
    end
end

function gather_mixed_cells(M::MCI, wd::WalkData, max_ind=0::Int)

    result = Set{MixedCell}()
    for c in keys(wd.walls)
        for (m, act_index, sgn) in wd.walls[c]
            !iszero(max_ind) && any(S -> any(i -> i > max_ind, S), m.inds) && continue
            push!(result, MixedCell(m.inds, M))
        end
    end

    return collect(result)
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

function id_matrix(n)
    return [i == j ? 1 : 0 for i in 1:n, j in 1:n]
end

function get_A_disc_equations(A::Matrix{Int})

    d, n = size(A)
    A_lift = vcat(A, id_matrix(n))
    A_lift = A_lift[1:size(A_lift, 1) - 1, :]
    edges = subsets([A_lift[:, i] for i in 1:n], 2)
    R, vars = polynomial_ring(QQ, vcat(["x$i" for i in 1:d], ["z$i" for i in 1:(n-1)]))
    x = vars[1:d]
    z = vars[d+1:end]
    s = sum(prod(map((i,j) -> i^j, vars, A_lift[:,k])) for k in 1:n)
    fs = [s; [v*derivative(s, v) for v in x]]
    return fs
end
