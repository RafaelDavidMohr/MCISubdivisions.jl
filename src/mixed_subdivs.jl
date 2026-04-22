# --- Functions for mixed subdivisions and homotopies --- #

function walk_homotopy!(w::WalkData)

    p0, p1 = w.p0, w.p1
    if p0.r == p1.r && p0.eps == p1.eps
        @info "no deformation, nothing to do"
        return
    end

    while !isempty(w.cells)
        mtbl = extract_min!(w.cells)
        if is_finished(mtbl)
            push!(w.finished_cells, mtbl.cell.inds)
            continue
        end
        new_cell_tables = mixed_cell_flip!(mtbl, w)
        for nc in new_cell_tables
            push!(w.cells, nc)
        end
    end
end

function deform_subdivision(A::Matrix{Int}, V::Matrix{C},
                            ms::Vector{MixedCellInds}, p0::DualVector,
                            p1::DualVector,
                            excluded_inds::Vector{Int}=Int[]) where C

    M = MCI(V, A)
    wd = WalkData(M, [MixedCell(m, M) for m in ms], p0, p1)
    walk_homotopy!(wd)
    return gather_mixed_cells(wd, excluded_inds)
end

function starting_system(A::Matrix{Int}, V::Matrix{C}, d::Vector{T}) where {C, T}
    n = size(A, 1)
    A_size = size(A, 2)
    no_multiset = allunique(i -> A[:, i], 1:A_size)

    # extended MCI
    p_start = rand(-LSIZE:LSIZE, A_size)
    col_inds = select_max_weight_columns(A, p_start)
    sd = @inbounds subdivision_of_points(transpose(A[:, col_inds]), -p_start[col_inds])
    A_start = copy(A)
    V_start = if C <: FqFieldElem
        F = parent(first(V))
        rand_arr_ff(F, n, A_size)
    else
        rand(-1000:1000, n, A_size)
    end
    init_mc = @inbounds [[col_inds[c]] for c in maximal_cells(sd)]

    return homotopy(A, V, A_start, V_start, init_mc, T.(p_start), d)
end

# works only if conv(A_target) ⊆ conv(A_start)
function homotopy(A_target::Matrix{Int}, V_target::Matrix{C},
                  A_start::Matrix{Int}, V_start::Matrix{C},
                  ms_start::Vector{MixedCellInds},
                  p_start::Vector{<:Integer},
                  p_target::Vector{<:Integer}) where C

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

function mixed_cell_flip!(mtbl::CellTable,
                          wd::WalkData)

    new_cell_tables = CellTable[]
    
    m = mtbl.cell
    @inbounds c = mtbl.walls[mtbl.i_min]
    ai = c.act_index
    ei = c.exchange_index
    
    # indices from which new mixed cell component can come
    new_inds = is_exchange(c) ? [ei] : @inbounds m.inds[ai + 1]

    M = wd.M
    
    Ss_new = Tuple{Int, Vector{Int}}[]
    @inbounds for (k, i) in enumerate(m.inds[ai])
        iszero(c.cfs[i]) && continue
        S_new = exchange(M.V, m.inds, ai, k, new_inds)
        if partial_sum_sign(S_new, c)
            push!(Ss_new, (k, S_new))
        end
    end

    inds = vcat(m.inds[ai], new_inds)
    for (k, S_new) in Ss_new
        sort!(S_new)
        S_new_next = sort(setdiff(inds, S_new))
        next_ind = is_exchange(c) ? ai + 1 : ai + 2
        is_split = length(S_new_next) > 1
        is_simple_exchange = is_exchange(c) && !is_split

        new_m_inds = if is_split
            @inbounds [m.inds[1:ai-1]..., S_new, S_new_next, m.inds[next_ind:end]...]
        else
            @inbounds [m.inds[1:ai-1]..., S_new, m.inds[next_ind:end]...]
        end

        MixedCell(new_m_inds, Vector{Int}[]) in wd.cells.seen && continue # this trick should work because of how we are hashing
        !is_affine_independent(M.A, new_m_inds) && continue

        new_loc_inds = if is_simple_exchange
            @inbounds si = m.inds[ai][k]
            @inbounds l = findfirst(ind -> ind == ei, m.loc_inds[ai])
            @inbounds new_loc_inds_ai = copy(m.loc_inds[ai])
            @inbounds new_loc_inds_ai[l] = si
            sort!(new_loc_inds_ai)
            @inbounds [i == ai ? new_loc_inds_ai : m.loc_inds[i] for i in 1:length(m.inds)]
        else
            known = @inbounds Dict([(i, m.loc_inds[i]) for i in 1:ai-1])
            shft = is_exchange(c) ? 0 : 1
            shft = is_split ? shft + 1 : shft
            for i in next_ind:length(m.inds)
                known[i + shft] = @inbounds m.loc_inds[i]
            end
            find_loc_indices!(M, new_m_inds, known)
        end

        @assert length(new_loc_inds) == length(new_m_inds)
        new_m = MixedCell(new_m_inds, new_loc_inds)

        new_tbl = if is_simple_exchange
            walls, i_min, t_min = new_cell_simple_exchange(mtbl, k, new_m)
            CellTable(new_m, walls, i_min, t_min)
        else
            t_last = cross_val_from_dp(c.dot0, c.dot1)
            walls, i_min, t_min = compute_active_walls!(new_m, M, wd.p0, wd.p1,
                                                        t_last)
            CellTable(new_m, walls, i_min, t_min)
        end
        push!(new_cell_tables, new_tbl)
    end
    return new_cell_tables
end

function exchange(V::Matrix{C}, m::MixedCellInds,
                  act_index::Int, k::Int,
                  new_inds::Vector{Int}) where C

    @inbounds candidate_indices = vcat(m[act_index][1:k-1], m[act_index][k+1:end],
                                       new_inds)
    @inbounds V_trunc_inds = vcat(candidate_indices,
                                  [idx[2:end] for idx in m[1:act_index-1]]...)
    @inbounds V_trunc = V[:, V_trunc_inds]
    cl = length(candidate_indices)
    return @inbounds candidate_indices[krnel_nz_entries(V_trunc, cl)]
end

function new_cell_simple_exchange(mtbl::CellTable{T},
                                  k::Int, new_m::MixedCell) where T

    m = mtbl.cell
    c = mtbl.walls[mtbl.i_min]
    t_last = mtbl.t_min
    
    ai = c.act_index
    si = m.inds[ai][k]

    walls = Hyperplane{T}[]

    dc0, dc1 = c.dot0, c.dot1
    t_min = one(DualNumber{T})
    i_min = 0
    @inbounds for (i, co) in enumerate(mtbl.walls)
        l1, l2 = if iszero(co.cfs[si]) || i == mtbl.i_min
            1, 0
        else
            l = lcm(c.cfs[si], co.cfs[si])
            div(l, co.cfs[si]), div(l, c.cfs[si]) # should always be ≠ 0
        end
        cfs_new = l1*co.cfs - l2*c.cfs
        dco0, dco1 = co.dot0, co.dot1
        d0_new = l1*dco0 - l2*dc0
        d1_new = l1*dco1 - l2*dc1
        new_ei = i == mtbl.i_min ? si : co.exchange_index
        t_min, i_min = add_hyperplane!(walls, new_m.inds, new_ei,
                                       co.act_index, t_min, i_min, t_last,
                                       cfs_new, d0_new, d1_new)
    end
    return walls, i_min, t_min
end

function compute_active_walls!(m::MixedCell,
                               M::MCI,
                               p0::DualVector{T},
                               p1::DualVector{T},
                               t_last::DualNumber{T}) where T

    k = length(m)
    n = ambient_dim(M)

    cayley_config = Matrix{Float64}(undef, n + k, 0)

    @inbounds for i in 1:k
        apd = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        for j in m.inds[i]
            cayley_config = hcat(cayley_config,
                                 vcat(M.A[:, j], apd))
        end
    end

    cayley_config_fac = factorize(cayley_config)
    d = Int(round(det(cayley_config_fac)))
    all_m_inds = vcat(m.inds...)

    walls = Hyperplane{T}[]
    t_min = one(DualNumber{T})
    i_min = 0
    @inbounds for i in 1:k
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
            end
            c_cfs_s = sparse(c_cfs)
            d0c, d1c = dot(p0, c_cfs_s), dot(p1, c_cfs_s)
            t_min, i_min = add_hyperplane!(walls, m.inds, ei, i,
                                           t_min, i_min,
                                           t_last, c_cfs_s, d0c, d1c)
        end
    end
    return walls, i_min, t_min
end

function add_hyperplane!(walls::Vector{Hyperplane{T}},
                         m::MixedCellInds,
                         ei::Int,
                         ai::Int,
                         t_min::DualNumber{T},
                         i_min::Int,
                         t_last::DualNumber{T},
                         cfs::SparseVector{Int, Int},
                         dc0::DualNumber,
                         dc1::DualNumber) where T

    @inbounds sgn = signbit(partial_sum(m[ai], cfs))
    c_new = Hyperplane{T}(cfs, dc0, dc1, ai, sgn, ei)
    t_c = cross_val_from_dp(c_new.dot0, c_new.dot1)
    push!(walls, c_new)
    if t_c <= t_last || one(DualNumber{T}) < t_c || t_c > t_min
        return t_min, i_min
    else
        return t_c, length(walls)
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
