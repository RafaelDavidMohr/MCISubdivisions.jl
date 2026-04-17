# --- Functions for mixed subdivisions and homotopies --- #

function walk_homotopy!(w::WalkData;
                        rr_counter=empty_rr_counter())

    trr = rr_counter.target_rr_count

    if trr >= 0
        ms = gather_mixed_cells(w)
        for m in ms
            add_mixed_cell!(rr_counter, m)
        end
        if rr_count(rr_counter) == trr
            @info "target real root count $(rr_counter.target_rr_count) reached"
            return true, zero(DualNumber)
        end
    end

    @info "starting homotopy"
    p0, p1 = w.p0, w.p1
    if p0.r == p1.r && p0.eps == p1.eps
        @info "no deformation, nothing to do"
        return true, zero(DualNumber)
    end

    while true
        c_int, t_cross = first_intersection!(w)
        if isnothing(c_int)
            @info "no intersection left, finished"
            return false, zero(DualNumber)
        end
        @info "intersection found"
        @info "crossing at $(t_cross)"
        walk_wall!(w, c_int, rr_counter)
        if trr >= 0 && rr_count(rr_counter) == trr
            @info "target real root count $(rr_counter.target_rr_count) reached"
            return true, t_cross
        end
        c_prev = c_int
    end
    return false, zero(DualNumber)
end

function deform_subdivision(A::Matrix{Int}, V::Matrix{C},
                            ms::Vector{MixedCellInds}, p0::DualVector,
                            p1::DualVector, target_rr_count::Int=-1) where C

    rr_counter = RRCounter(A, V, target_rr_count)

    M = MCI(V, A)
    wd = WalkData(M, [MixedCell(m, M) for m in ms], p0, p1)
    cnt_reached, t_cross = walk_homotopy!(wd, rr_counter = rr_counter)
    return cnt_reached, t_cross, gather_mixed_cells(wd)
end

function starting_system(A::Matrix{Int}, V::Matrix{FqFieldElem}, d::Vector{Int})
    n = size(A, 1)
    A_size = size(A, 2)
    no_multiset = allunique(i -> A[:, i], 1:A_size)

    # extended MCI
    F = parent(first(V))
    p_start = rand(-10000:10000, A_size)
    @info "Building start system"
    col_inds = select_max_weight_columns(A, p_start)
    sd = subdivision_of_points(transpose(A[:, col_inds]), -p_start[col_inds])
    A_start = copy(A)
    V_start = rand_arr_ff(F, n, A_size)
    init_mc = [[col_inds[c]] for c in maximal_cells(sd)]

    return homotopy(A, V, A_start, V_start, init_mc, p_start, d)
end

# works only if conv(A_target) ⊆ conv(A_start)
function homotopy(A_target::Matrix{Int}, V_target::Matrix{FqFieldElem},
                  A_start::Matrix{Int}, V_start::Matrix{FqFieldElem},
                  ms_start::Vector{MixedCellInds},
                  p_start::Vector{Int},
                  p_target::Vector{Int})

    A_ext = hcat(A_target, A_start)
    V_ext = hcat(V_target, V_start)

    M = MCI(V_ext, A_ext)

    A_size = size(A_target, 2)

    # set up path
    A_ext_size = size(A_ext, 2)
    p0 = forgetful_lift(A_ext_size, collect(1:A_size), p_start)
    p1 = forgetful_lift(A_ext_size, collect(A_size+1:A_ext_size), p_target)

    # initial mixed cells
    init_mc = [MixedCell([ci .+ A_size for ci in c], M) for c in ms_start]

    wd = WalkData(M, init_mc, p0, p1)

    return wd
end

# --- Functions related to mixed cell cones --- #

function walk_wall!(wd::WalkData, c::Hyperplane, rr_counter::RRCounter)
    active_mc_data = wd.walls[c]
    delete!(wd.walls, c) 
    new_ms_inds = Set{MixedCellInds}()
    new_ms = MixedCell[]
    @info "$(length(active_mc_data)) cells to flip"
    for (m, act_index, sgn) in active_mc_data
        m_circuits = delete_mixed_cell!(wd, m, rr_counter)
        mixed_cell_flip!(m, c, wd, act_index, sgn, m_circuits, new_ms_inds, new_ms)
    end
    @info "$(length(new_ms_inds)) new mixed cells"
    if rr_counter.target_rr_count >= 0
        for new_mc in new_ms
            add_mixed_cell!(rr_counter, new_mc.inds)
        end
    end
end

function mixed_cell_flip!(m::MixedCell, c::Hyperplane, wd::WalkData,
                          act_index::Int,
                          sgn::Bool,
                          m_circuits::Vector{Tuple{Hyperplane, Int}},
                          new_ms_inds::Set{MixedCellInds},
                          new_ms::Vector{MixedCell})

    # indices from which new mixed cell component can come
    new_inds = setdiff(nz_inds(c), vcat(m.inds...))
    is_exchange = !isempty(new_inds)
    if is_exchange
        @assert isone(length(new_inds))
    else
        new_inds = m.inds[act_index+1]
    end

    M = wd.M
    
    Ss_new = Tuple{Int, Vector{Int}}[]
    F = prime_field_V(M)
    for (k, i) in enumerate(m.inds[act_index])
        iszero(c.cfs[i]) && continue # correct?
        S_new = exchange(M.V, m.inds, act_index, k, new_inds)
        if partial_sum_sign(S_new, c, sgn)
            push!(Ss_new, (k, S_new))
        end
    end

    inds = vcat(m.inds[act_index], sort(new_inds))
    for (k, S_new) in Ss_new
        sort!(S_new)
        S_new_next = sort(setdiff(inds, S_new))
        next_ind = is_exchange ? act_index + 1 : act_index + 2
        new_m_inds = if length(S_new_next) > 1
            [m.inds[1:act_index-1]..., S_new, S_new_next,
             m.inds[next_ind:end]...]
        else
            [m.inds[1:act_index-1]..., S_new,
             m.inds[next_ind:end]...]
        end
        new_m_inds in new_ms_inds && continue
        if is_affine_independent(M.A, new_m_inds)
            is_simple_exchange = is_exchange && length(S_new_next) <= 1
            new_mc = if is_simple_exchange
                new_ind = first(new_inds)
                new_cell_simple_exchange!(m, act_index, k, new_ind,
                                          new_m_inds, wd, c, m_circuits)
            else
                new_mc = MixedCell(new_m_inds, M)
                compute_active_walls!(new_mc, M, wd.cells, wd.walls, wd.p0, wd.p1, wd.t_curr)
                new_mc
            end
            push!(new_ms, new_mc)
            push!(new_ms_inds, new_m_inds)
        end
    end
end

function exchange(V::Matrix{C}, m::MixedCellInds, act_index::Int, k::Int, new_inds::Vector{Int}) where C
    candidate_indices = vcat(m[act_index][1:k-1], m[act_index][k+1:end],
                             new_inds)
    V_trunc_inds = vcat(candidate_indices,
                        [idx[2:end] for idx in m[1:act_index-1]]...)
    V_trunc = V[:, V_trunc_inds]
    F = parent(first(V))
    K = kernel(matrix(F, V_trunc), side = :right)
    @assert isone(size(K, 2))
    cl = length(candidate_indices)
    return candidate_indices[findall(j -> !iszero(K[j, :]), 1:cl)]
end

function new_cell_simple_exchange!(m::MixedCell, act_index::Int, k::Int, # position in mixed cell
                                   new_ind::Int,
                                   new_m_inds::Vector{Vector{Int}},
                                   wd::WalkData, c::Hyperplane,
                                   m_circuits::Vector{Tuple{Hyperplane, Int}})

    ei = m.inds[act_index][k]
    new_loc_inds = sort(vcat(setdiff(m.loc_inds[act_index], [new_ind]), [ei]))
    
    new_m = MixedCell(new_m_inds, [i == act_index ? new_loc_inds : m.loc_inds[i] for i in 1:length(m.inds)])

    wd.cells[new_m] = Hyperplane[]
    c_cfs = c.cfs
    for (co, ai) in m_circuits
        co == c && continue
        new_c = if iszero(co.cfs[ei])
            co
        else
            co_cfs = co.cfs
            l = lcm(c_cfs[ei], co_cfs[ei])
            l1 = div(l, co_cfs[ei]) 
            l2 = div(l, c_cfs[ei]) # should always be ≠ 0
            Hyperplane(l1*co_cfs - l2*c_cfs, l1*co.dot0 - l2*c.dot0, l1*co.dot1 - l2*c.dot1)
        end
        add_hyperplane!(new_m, new_c, ai, wd.cells, wd.walls, wd.t_curr)
    end

    return new_m
end

function compute_active_walls!(m::MixedCell,
                               M::MCI,
                               cells::Dict{MixedCell, Vector{Hyperplane}},
                               walls::Dict{Hyperplane, Vector{Tuple{MixedCell, Int, Bool}}},
                               p0::DualVector,
                               p1::DualVector,
                               t_cross::Union{Nothing, DualNumber})

    k = length(m)
    n = ambient_dim(M)

    cayley_config = Matrix{Float64}(undef, n + k, 0)

    for i in 1:k
        apd = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        for j in m.inds[i]
            cayley_config = hcat(cayley_config,
                                 vcat(M.A[:, j], apd))
        end
    end

    cayley_config_fac = factorize(cayley_config)
    d = Int(round(det(cayley_config_fac)))
    all_m_inds = vcat(m.inds...)

    cells[m] = Hyperplane[]
    for i in 1:k
        apd = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        for j in cayley_indices(m, i)
            c_cfs = zeros(Int, size(M.A, 2))
            idx = findfirst(l -> M.A[:, l] == M.A[:, j], m.inds[i])
            if !isnothing(idx)
                c_cfs[m.inds[i][idx]] = 1
                c_cfs[j] = -1
            else
                rs = vcat(M.A[:, j], apd)
                sol = cayley_config_fac \ rs
                for (l, ind) in enumerate(all_m_inds)
                    c_cfs[ind] = Int(round(d * sol[l]))
                end
                c_cfs[j] -= d
                g = gcd(c_cfs)
                c_cfs = [div(cf, g) for cf in c_cfs]
            end
            c = Hyperplane(c_cfs, p0, p1)
            add_hyperplane!(m, c, i, cells, walls, t_cross)
        end
    end
end

function add_hyperplane!(m::MixedCell, c::Hyperplane,
                         i::Int,
                         cells::Dict{MixedCell, Vector{Hyperplane}},
                         walls::Dict{Hyperplane, Vector{Tuple{MixedCell, Int, Bool}}},
                         t_cross::Union{Nothing, DualNumber})

    t_c = c.cross_val
    (lt_dual(t_c, zero(DualNumber)) || lt_dual(one(DualNumber), t_c)) && return
    if !isnothing(t_cross) && (t_c == t_cross || lt_dual(t_c, t_cross))
        return
    end
    push!(cells[m], c)
    sgn = signbit(partial_sum(m.inds[i], c))
    if haskey(walls, c)
        if all(mtp -> mtp[1] != m, walls[c])
            push!(walls[c], (m, i, sgn))
        end
    else
        walls[c] = [(m, i, sgn)]
    end
end

# --- Mixed cell data --- #

function outer_normal_vector(A::Matrix{Int}, m::MixedCellInds,
                             d::Vector{C}) where C

    A_lifted = vcat(A, transpose(d))
    K = normal_space(A_lifted, m)
    @assert isone(size(K, 2)) "mixed cell does not lift to a hyperplane"
    v = K[:, 1]
    v *= v[end]^(-1)
    n = length(v)
    return v[1:n-1]
end

function outer_normal_vector(M::MCI, m::MixedCell,
                             d::Vector{C}) where C

    return outer_normal_vector(M.A, m, d)
end
