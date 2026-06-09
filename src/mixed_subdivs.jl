# --- Functions for mixed subdivisions and homotopies --- #

function walk_homotopy!(w::WalkData)

    p0, p1 = w.p0, w.p1
    if [poi.r for poi in p0] == [p1i.r for p1i in p1] # TEMPORARY
        @info "no deformation, nothing to do"
        while !isempty(w.cells)
            mtbl = extract_min!(w.cells)
            push!(w.finished_cells, mtbl.cell.inds)
        end
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

function starting_system(A::Matrix{Int}, V::Matrix{C}, d::Vector{Int}) where C
    n = size(A, 1)
    A_size = size(A, 2)

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

    return homotopy(A, V, A_start, V_start, init_mc, p_start, d)
end

# works only if conv(A_target) ⊆ conv(A_start)
function homotopy(A_target::Matrix{Int}, V_target::Matrix{C},
                  A_start::Matrix{Int}, V_start::Matrix{C},
                  ms_start::Vector{MixedCellInds},
                  p_start::Vector{Int},
                  p_target::Vector{Int}) where C

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
        iszero(c.cfs_modP[i]) && continue
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
            t_last = mtbl.t_min
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

function new_cell_simple_exchange(mtbl::CellTable,
                                  k::Int, new_m::MixedCell)

    m = mtbl.cell
    c = mtbl.walls[mtbl.i_min]
    t_last = mtbl.t_min
    
    ai = c.act_index
    si = m.inds[ai][k]

    walls = Hyperplane[]

    dc0, dc1 = c.dot0, c.dot1
    dc0P, dc1P = c.dot0_modP, c.dot1_modP
    t_min = one(OrderDual)
    i_min = 0
    @inbounds for (i, co) in enumerate(mtbl.walls)
        mult, multP = if iszero(co.cfs_modP[si]) || i == mtbl.i_min
            0.0, zero(FPNum)
        else
            c.cfs[si]^(-1) * co.cfs[si], inv(c.cfs_modP[si]) * co.cfs_modP[si]
        end
        cfs_new = co.cfs - mult*c.cfs
        cfs_new[si] = 0.0
        cfs_newP = co.cfs_modP - multP*c.cfs_modP
        dco0, dco1 = co.dot0, co.dot1
        d0_new = dco0 - mult*dc0
        d1_new = dco1 - mult*dc1
        dco0P, dco1P = co.dot0_modP, co.dot1_modP
        d0_newP = dco0P - multP*dc0P
        d1_newP = dco1P - multP*dc1P
        new_ei = i == mtbl.i_min ? si : co.exchange_index
        t_min, i_min = add_hyperplane!(walls, new_m.inds, new_ei,
                                       co.act_index, t_min, i_min, t_last,
                                       cfs_new, cfs_newP,
                                       d0_new, d1_new,
                                       d0_newP, d1_newP)
    end
    return walls, i_min, t_min
end

function compute_active_walls!(m::MixedCell,
                               M::MCI,
                               p0::DualVector{Int},
                               p1::DualVector{Int},
                               t_last::OrderDual)

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

    walls = Hyperplane[]
    t_min = one(OrderDual)
    i_min = 0
    @inbounds for i in 1:k
        apd = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        for j in cayley_indices(m, i)
            ei = i == k || j ∉ m.inds[i+1] ? j : 0
            c_cfs = zeros(Float64, size(M.A, 2))
            c_cfs_modP = zeros(FPNum, size(M.A, 2))
            idx = findfirst(l -> M.A[:, l] == M.A[:, j], m.inds[i])
            if !isnothing(idx)
                c_cfs[m.inds[i][idx]] = 1.0
                c_cfs[j] = -1.0
                c_cfs_modP[m.inds[i][idx]] = one(FPNum)
                c_cfs_modP[j] = -one(FPNum)
            else
                rs = vcat(M.A[:, j], apd)
                sol = cayley_config_fac \ rs
                for (l, ind) in enumerate(all_m_inds)
                    c_cfs[ind] = sol[l]
                    c_cfs_modP[ind] = FPNum(mod(Int(round(d * sol[l])), PRIME))
                end
                c_cfs[j] -= 1
                c_cfs_modP[j] -= FPNum(d % PRIME)
            end
            c_cfs_s = sparse(c_cfs)
            c_cfs_modP_s = sparse(c_cfs_modP)
            d0c, d1c = dot(p0, c_cfs_s), dot(p1, c_cfs_s)
            d0c_modP, d1c_modP = dot(p0, c_cfs_modP_s), dot(p1, c_cfs_modP_s)
            t_min, i_min = add_hyperplane!(walls, m.inds, ei, i,
                                           t_min, i_min,
                                           t_last, c_cfs_s, c_cfs_modP_s,
                                           d0c, d1c,
                                           d0c_modP, d1c_modP)
        end
    end
    return walls, i_min, t_min
end

function add_hyperplane!(walls::Vector{Hyperplane},
                         m::MixedCellInds,
                         ei::Int,
                         ai::Int,
                         t_min::OrderDual,
                         i_min::Int,
                         t_last::OrderDual,
                         cfs::SparseVector{Float64, Int},
                         cfs_modP::SparseVector{FPNum, Int},
                         dc0::DualNumber,
                         dc1::DualNumber,
                         dc0_modP::DualNumber,
                         dc1_modP::DualNumber)

    @inbounds sgn = signbit(partial_sum(m[ai], cfs))
    c_new = Hyperplane(cfs, cfs_modP,
                       dc0, dc1,
                       dc0_modP, dc1_modP,
                       ai, sgn, ei)
    t_c = cross_val(c_new)
    push!(walls, c_new)
    if t_c <= t_last || one(OrderDual) < t_c || t_c > t_min
        return t_min, i_min
    else
        return t_c, length(walls)
    end
end

# --- Mixed cell data --- #

function primitive_normal_vector(A::Matrix{Int}, m::MixedCellInds)

    K = normal_space(A, m)
    @assert isone(size(K, 2)) "mixed cell does not lift to a hyperplane"
    v = K[:, 1]
    l = lcm(denominator.(v))
    res = Int.(l * v)
    g = gcd(res)
    res = res .÷ g
    test_ind = findfirst(idx -> all(mi -> idx ∉ mi, m), 1:size(A, 2))
    if dot(res, A[:, test_ind]) > dot(res, A[:, first(first(m))])
        return -res
    else
        return res
    end
end

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
