# --- Functions for mixed subdivisions and homotopies --- #

function walk_homotopy!(w::WalkData)

    @info "starting homotopy"
    p0, p1 = w.p0, w.p1
    if p0.r == p1.r && p0.eps == p1.eps
        @info "no deformation, nothing to do"
        return
    end

    while !isempty(w.cells)
        @info "$(length(w.cells)) cells left"
        mtbl = popfirst!(w.cells)
        if is_finished(mtbl)
            @info "cell finished"
            push!(w.finished_cells, mtbl.cell.inds)
            continue
        end
        new_cell_tables = mixed_cell_flip!(mtbl, w)
        append!(w.cells, new_cell_tables)
    end
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

function mixed_cell_flip!(mtbl::MixedCell,
                          wd::WalkData)

    new_cell_tables = CellTable[]
    
    m = mtbl.cell
    c = mtbl.walls[mtbl.i_min]
    ai = c.act_index
    
    # indices from which new mixed cell component can come
    new_inds = is_exchange(c) ? [c.exchange_index] : m.inds[ai + 1]

    M = wd.M
    
    Ss_new = Tuple{Int, Vector{Int}}[]
    F = prime_field_V(M)
    for (k, i) in enumerate(m.inds[ai])
        iszero(c.cfs[i]) && continue # correct?
        S_new = exchange(M.V, m.inds, act_index, k, new_inds)
        if partial_sum_sign(S_new)
            push!(Ss_new, (k, S_new))
        end
    end

    inds = vcat(m.inds[ai], new_inds)
    for (k, S_new) in Ss_new
        sort!(S_new)
        S_new_next = sort(setdiff(inds, S_new))
        next_ind = is_exchange(c) ? act_index + 1 : act_index + 2
        new_m_inds = if length(S_new_next) > 1
            [m.inds[1:act_index-1]..., S_new, S_new_next,
             m.inds[next_ind:end]...]
        else
            [m.inds[1:act_index-1]..., S_new,
             m.inds[next_ind:end]...]
        end
        if is_affine_independent(M.A, new_m_inds)
            is_simple_exchange = is_exchange(c) && length(S_new_next) <= 1
            new_tbl = if is_simple_exchange
                new_m, walls, i_min = new_cell_simple_exchange(mtbl, k, new_m_inds)
                CellTable(new_m, walls, i_min)
            else
                new_m = MixedCell(new_m_inds, M)
                t_last = cross_val_from_dp(c.dot0, c.dot1)
                walls, i_min = compute_active_walls!(new_mc, M, wd.p0, wd.p1, t_last)
                CellTable(new_m, walls, i_min)
            end
            push!(new_cell_tables, new_tbl)
        end
    end
    return new_cell_tables
end

function exchange(V::Matrix{C}, m::MixedCellInds,
                  act_index::Int, k::Int,
                  new_inds::Vector{Int}) where C
    
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

function new_cell_simple_exchange(mtbl::CellTable, k::Int, # position in mixed cell
                                  new_m_inds::Vector{Vector{Int}})

    m = mtbl.cell
    c = mtbl.walls[mtbl.i_min]
    t_last = cross_val_from_dp(c.dot0, c.dot1)
    
    si = m.inds[act_index][k]
    ei = c.exchange_index
    ai = c.act_index
    l = findfirst(ind -> ind == ei, m.loc_inds[ai])
    new_loc_inds = m.loc_inds[ai]
    new_loc_inds[l] = si
    sort!(new_loc_inds)
    
    new_m = MixedCell(new_m_inds, [i == ai ? new_loc_inds : m.loc_inds[i] for i in 1:length(m.inds)])

    walls = Hyperplane[]

    dc0, dc1 = c.dot0, c.dot1
    t_min = one(DualNumber)
    i_min = 0
    for (i, co) in enumerate(mtbl.walls)
        i == c_ind && continue # THIS MAY BREAK IF THERE IS A DUPLICATE
        if iszero(co.cfs[si])
            push!(walls, co)
            t_co = cross_val_from_dp(co.dot0, co.dot1)
            if t_co < t_min
                t_min = t_co
                i_min = length(walls)
            end
        else
            l = lcm(c.cfs[si], co.cfs[si])
            l1 = div(l, co.cfs[si]) 
            l2 = div(l, c.cfs[si]) # should always be ≠ 0
            cfs_new = l1*co.cfs - l2*c.cfs
            dco0, dco1 = co.dot0, co.dot1
            d0_new = l1*dco0 - l2*dc0
            d1_new = l1*dco1 - l2*dc1
            t_min, i_min = add_hyperplane!(walls, new_m.inds, co.exchange_index,
                                           co.act_index, t_min, i_min, t_last,
                                           cfs_new, d0_new, d1_new)
        end
    end
    return new_m, walls, i_min
end

function compute_active_walls!(m::MixedCell,
                               M::MCI,
                               p0::DualVector,
                               p1::DualVector,
                               t_last::DualNumber)

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

    walls = Hyperplane[]
    t_min = one(DualNumber)
    i_min = 0
    for i in 1:k
        apd = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        for j in cayley_indices(m, i)
            ei = i == k || j ∉ m.inds[i+1] ? j : 0
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
            d0c, d1c = dot(p0, cfs), dot(p1, cfs)
            t_min, i_min = add_hyperplane!(walls, m.inds, ei, i,
                                           t_min, i_min,
                                           t_last, c, d0c, dc1)
        end
    end
    return walls, i_min
end

function add_hyperplane!(walls::Vector{Hyperplane},
                         m::MixedCellInds,
                         ei::Int,
                         ai::Int,
                         t_min::DualNumber,
                         i_min::Int,
                         t_last::DualNumber,
                         cfs::SparseVector{Int, Int},
                         dc0::DualNumber,
                         dc1::DualNumber)

    t_c = cross_val_from_dp(dc0, dc1)
    if t_c <= t_last || one(DualNumber) < t_c
        return t_min, i_min
    end
    sgn = signbit(partial_sum(m[i], c_new))
    c_new = Hyperplane(cfs, dc0, dc1, ai, sgn, ei)
    push!(walls, c_new)
    if t_c < t_min
        return t_c, length(walls)
    else
        return t_min, i_min
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
