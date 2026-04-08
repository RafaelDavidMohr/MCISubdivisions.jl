# --- mixed cells --- #

function MixedCell(inds::Vector{Vector{Int}}, M::MCI)
    for ind in inds
        sort!(ind)
    end
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

# assume that V is echelonized
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
    return vol(m.inds, A)
end

function vol(m::Vector{Vector{Int}}, A::Matrix{Int})
    mat = hcat([linear_span(A, S) for S in m]...)
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
        return sort(setdiff(m.loc_inds[j], m.inds[j]))
    end
    res = setdiff(m.loc_inds[j], vcat(m.inds[j], m.loc_inds[j+1]))
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

function linear_span(A::Matrix{C}, inds::Vector{Int}) where C
    a0 = A[:, first(inds)]
    L = Matrix{C}(undef, size(A, 1), length(inds) - 1)
    for (j, i) in enumerate(inds[2:end])
        L[:, j] = A[:, i] - a0
    end
    return L
end

function is_affine_independent(A::Matrix{Int}, inds::Vector{Vector{Int}})
    mat = hcat([linear_span(A, S) for S in inds]...)
    return !iszero(round(det(mat))) 
end

function normal_space(A::Matrix{C}, m::Vector{Vector{Int}}) where C

    n = size(A, 1)
    eqns = Matrix{C}(undef, n, 0)
    for S in m
        a = A[:, first(S)]
        for i in S[2:end]
            eqns = hcat(eqns, A[:, i] - a)
        end
    end

    if C <: Float64
        return nullspace(eqns)
    end

    F = if C <: Int
        QQ
    else
        parent(first(A))
    end

    return permutedims(Matrix(kernel(matrix(F, eqns))))
end

# --- circuits --- #

function nz_inds(c::Hyperplane)
    return c.nzinds
end

function partial_sum(inds::Vector{Int}, c::Hyperplane)
    return sum(c.cfs[inds])
end

function partial_sum_sign(inds::Vector{Int}, c::Hyperplane, sgn::Bool)
    res = partial_sum(inds, c)
    return res != 0 && (sgn ? res > 0 : res < 0) # TODO: check if 0 is allowed here
end

function LinearAlgebra.dot(v::Vector{Int}, c::Hyperplane)
    res = 0
    for i in nz_inds(c)
        res += v[i] * c.cfs[i]
    end
    return res
end

function LinearAlgebra.dot(v::Vector{Float64}, c::Hyperplane)
    res_fl = 0.0
    for i in nz_inds(c)
        res_fl += v[i] * c.cfs[i]
    end
    return res_fl
end

function LinearAlgebra.dot(l::DualVector, c::Hyperplane)
    return DualNumber(round(dot(l.r, c)), round(dot(l.eps, c)))
end

# --- homotopy paths --- #
    
function first_intersection!(l0::DualVector, l1::DualVector,
                             hyperplanes::AbstractSet{Hyperplane},
                             last_h::Union{Nothing, Hyperplane},
                             wd::WalkData)

    last_t = if isnothing(last_h)
        zero(DualNumber)
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
        if lt_dual(t_c, last_t) || lt_dual(one(DualNumber), t_c)
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

    if iszero(l0d.r) && iszero(l1d.r)
        return DualNumber(l0d.eps * denom.eps^(-1), 0.0)
    elseif iszero(denom.r)
        return DualNumber(2.0, 0.0)
    else
        return l0d * inv(denom)
    end
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
    return DualVector(d)
end

function forgetful_lift(A_size::Int, forget_inds::Vector{Int}, d_eps::Vector{Int})
    d = Vector{Int}(undef, A_size)
    d_eps_new = Vector{Int}(undef, A_size)
    cnt = 1
    for i in 1:A_size
        if i in forget_inds
            d[i] = 0
            d_eps_new[i] = rand(-10000:10000)
        else
            d[i] = 1
            d_eps_new[i] = d_eps[cnt]
            cnt += 1
        end
    end
    return DualVector(d, d_eps_new)
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

function delete_mixed_cell!(wd::WalkData, m::MixedCell, rr_counter::RRCounter)
    for c in keys(wd.walls)
        filter!(((m0, i, j),) -> m != m0, wd.walls[c])
        if isempty(wd.walls[c])
            delete!(wd.walls, c)
        end
    end
    delete!(rr_counter.cnt, m)
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

function id_matrix(n)
    return [i == j ? 1 : 0 for i in 1:n, j in 1:n]
end

function get_A_disc_equations(A::Matrix{Int}, dehom=1)

    d, n = size(A)
    A_lift = vcat(A, id_matrix(n))
    A_lift = A_lift[1:size(A_lift, 1) .!= (dehom + d), :]
    R, vars = polynomial_ring(QQ, vcat(["x$i" for i in 1:d], ["z$i" for i in 1:(n-1)]))
    x = vars[1:d]
    z = vars[d+1:end]
    s = sum(rand(-1000:1000)*prod(map((i,j) -> i^j, vars, A_lift[:,k])) for k in 1:n)
    fs = [s; [v*derivative(s, v) for v in x]]
    return fs
end
