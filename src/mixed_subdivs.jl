# --- Functions for mixed subdivisions and homotopies --- #

function walk_homotopy!(w::WalkData, p0::DualVector, p1::DualVector;
                        rr_counter=empty_rr_counter())

    trr = rr_counter.target_rr_count

    if trr >= 0
        ms = gather_mixed_cells(w.M, w)
        for m in ms
            add_mixed_cell!(rr_counter, m)
        end
        if rr_count(rr_counter) == trr
            @info "target real root count $(rr_counter.target_rr_count) reached"
            return true, zero(DualNumber)
        end
    end

    @info "starting homotopy"
    if p0.r == p1.r && p0.eps == p1.eps
        @info "no deformation, nothing to do"
        return true, zero(DualNumber)
    end

    c_prev = nothing
    while true
        c_int, t_cross = first_intersection!(p0, p1, keys(w.walls), c_prev, w)
        if isnothing(c_int)
            @info "no intersection left, finished"
            return false, zero(DualNumber)
        end
        @info "intersection found"
        @info "crossing at $(t_cross)"
        @info "$(length(w.walls[c_int])) mixed cells to flip"
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
                            ms::Vector{MixedCell}, p0::DualVector,
                            p1::DualVector, target_rr_count::Int) where C

    rr_counter = RRCounter(A, V, target_rr_count)

    M = MCI(V, A)
    wd = WalkData(M, ms)
    cnt_reached, t_cross = walk_homotopy!(wd, p0, p1, rr_counter = rr_counter)
    return cnt_reached, t_cross, gather_mixed_cells(M, wd)
end

function starting_system(A::Matrix{Int}, V::Matrix{FqFieldElem})
    n = size(A, 1)
    A_size = size(A, 2)
    no_multiset = allunique(i -> A[:, i], 1:A_size)

    # extended MCI
    A_ext = copy(A)
    V_ext = copy(V)
    F = parent(first(V))
    if no_multiset
        @info "Choosing polyhedral starting system"
        A_ext = hcat(A_ext, A)
        V_ext = hcat(V_ext, rand_arr_ff(F, n, A_size))
    else
        @info "Choosing total degree starting system"
        max_deg = maximum(i -> sum(A[:, i]), 1:A_size)
        for i in 0:n
            A_ext = hcat(A_ext, [j == i ? max_deg : 0 for j in 1:n])
            V_ext = hcat(V_ext, rand_arr_ff(F, n))
        end
    end
    M = MCI(V_ext, A_ext)

    # set up path
    A_ext_size = size(A_ext, 2)
    p0 = forgetful_lift(A_ext_size, collect(1:A_size))
    p1 = forgetful_lift(A_ext_size, collect(A_size+1:A_ext_size))

    # initial mixed cells
    init_mc = if no_multiset
        sd = subdivision_of_points(transpose(A), -p0.eps[A_size+1:end])
        [MixedCell([c .+ A_size], M) for c in maximal_cells(sd)]
    else
        [MixedCell([collect(A_size+1:A_ext_size)], M)]
    end

    wd = WalkData(M, init_mc)

    return A_ext, wd, p0, p1
end

# --- Functions related to mixed cell cones --- #

function walk_wall!(wd::WalkData, c::Hyperplane, rr_counter::RRCounter)
    active_mc_data = wd.walls[c]
    delete!(wd.walls, c) 
    cnt = 0
    for (m, act_index, sgn) in active_mc_data
        delete_mixed_cell!(wd, m, rr_counter)
        new_mixed_cells = mixed_cell_flip(m, c, wd.M, act_index, sgn)
        cnt += length(new_mixed_cells)
        for new_mc in new_mixed_cells
            compute_active_walls!(new_mc, wd.M, wd.walls)
            if rr_counter.target_rr_count >= 0
                add_mixed_cell!(rr_counter, new_mc)
            end
        end
    end
    @info "$(cnt) new mixed cells (possible duplicates)"
end

function mixed_cell_flip(m::MixedCell, c::Hyperplane, M::MCI, act_index::Int, sgn::Bool)

    # indices from which new mixed cell component can come
    # todo: if this doesnt work check if this is correct
    new_inds = setdiff(nz_inds(c), vcat(m.inds...))
    is_exchange = !isempty(new_inds)
    is_exchange && @assert isone(length(new_inds))

    Ss_new = Vector{Int}[]
    for (k, i) in enumerate(m.inds[act_index])
        (iszero(c.cfs[i]) || signbit(c.cfs[i]) != sgn) && continue
        candidate_indices = vcat(m.inds[act_index][1:k-1], m.inds[act_index][k+1:end],
                                 new_inds)
        V_trunc_inds = vcat(candidate_indices, m.inds[1:act_index-1]...)
        V_trunc = M.V[:, V_trunc_inds]
        cl = length(candidate_indices)
        K_proj = project_kernel(V_trunc, collect(1:cl))
        @assert isone(size(K_proj, 2))
        S_new = candidate_indices[findall(j -> !iszero(K_proj[j, :]), 1:cl)]
        push!(Ss_new, S_new)
    end

    inds = is_exchange ? vcat(m.inds[act_index], sort(new_inds)) : vcat(m.inds[act_index], m.inds[act_index+1])
    new_mixed_cells = MixedCell[]
    for S_new in Ss_new
        S_new_next = sort(setdiff(inds, S_new))
        next_ind = is_exchange ? act_index + 1 : act_index + 2
        new_ms_inds = if length(S_new_next) > 1
            [m.inds[1:act_index-1]..., S_new, S_new_next,
             m.inds[next_ind:end]...]
        else
            [m.inds[1:act_index-1]..., S_new,
             m.inds[next_ind:end]...]
        end
        if is_affine_independent(M.A, new_ms_inds)
            push!(new_mixed_cells, MixedCell(new_ms_inds, M))
        end
    end

    return new_mixed_cells
end

function compute_active_walls!(m::MixedCell,
                               M::MCI,
                               walls::Dict{Hyperplane, Set{Tuple{MixedCell, Int, Bool}}})

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

    cayley_config = factorize(cayley_config)
    d = Int(round(det(cayley_config)))
    all_m_inds = vcat(m.inds...)

    for i in 1:k
        apd = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        for j in cayley_indices(m, i)
            rs = vcat(M.A[:, j], apd)
            sol = cayley_config \ rs
            c_cfs = zeros(Int, size(M.A, 2))
            for (l, ind) in enumerate(all_m_inds)
                c_cfs[ind] = Int(round(d * sol[l]))
            end
            c_cfs[j] -= d
            c = Hyperplane(c_cfs)
            sgn = signbit(first(partial_sum(m.inds[i], c)))
            add_to_dict!(walls, c, (m, i, sgn))
        end
    end
end

# --- Mixed cell checking/computation for testing --- #

function outer_normal_vector(A::Matrix{Int}, m::MixedCell,
                             d::Vector{C}) where C

    

    A_lifted = vcat(A, transpose(d))
    K = normal_space(A_lifted, m.inds)
    @assert isone(size(K, 2)) "mixed cell does not lift to a hyperplane"
    K *= K[end, 1]^(-1)
    return K[1:n, 1]
end

function outer_normal_vector(M::MCI, m::MixedCell,
                             d::Vector{C}) where C

    A_modP_lifted = vcat(M.A_modP, transpose((F).(d)))
    A_Fl_lifted = vcat(M.A_Fl, transpose((Float64).(d)))

    K_modP = normal_space(A_modP_lifted, m.inds)
    @assert isone(size(K_modP, 2)) "mixed cell does not lift to a hyperplane"
    K_modP *= K_modP[end, 1]^(-1)

    A_Fl_lifted = vcat(M.A_Fl, transpose((Float64).(d)))
    K_Fl = normal_space(A_Fl_lifted, m.inds)
    K_Fl *= K_Fl[end, 1]^(-1)

    return K_modP[1:n, 1], K_Fl[1, 1:n]
end

function find_dual_tropical_root(M::MCI, d::Vector{QQFieldElem},
                                 w::Vector{QQFieldElem})

    result = Vector{Int}[]

    A_card = size(M.A_modP, 2)
    A_indices = collect(1:A_card)
    n = ambient_dim(M)
    w_fl = (Float64).(w)
    d_fl = (Float64).(d)
    sort!(A_indices, by = i -> dot(w_fl, M.A_Fl[:, i]) + d_fl[i], rev = true)

    F = prime_field_A(M)
    w_modP = (F).(w)
    d_modP = (F).(d)

    prev_deg = dot(w_modP, M.A_modP[:, first(A_indices)]) + d_modP[first(A_indices)]
    Sj = Int[]
    result_codim = 0

    for (j, i) in enumerate(A_indices)
        result_codim == n && break
        new_deg = dot(w_modP, M.A_modP[:, i]) + d_modP[i] 
        if new_deg == prev_deg
            push!(Sj, i)
        end
        if j == length(A_indices) || new_deg != prev_deg
            result_codim += length(Sj) - 1
            push!(result, copy(Sj))
            Sj = [i]
        end
        prev_deg = new_deg
    end
            
    return MixedCell(result, M)
end

# checks if m ∪ {S} is a partial mixed cell
# function is_partial_mixed_cell(M::MCI, m::MixedCell)

#     M_curr = RelativeMCI(M)

#     for j in 1:length(m.inds)
#         !is_partial_mixed_cell(M_curr, m.inds[j]) && return false
#         M_curr = localize(M_curr, m.inds[j], m.loc_inds[j])
#     end

#     return true
# end

# function is_partial_mixed_cell(M::RelativeMCI,
#                                S::Vector{Int})

#     S_rel = M.base_to_rel[S]
#     A = M.A_rel
#     V = M.V_rel
#     if is_affine_independent(A, S_rel)
#         F = parent(first(V))
#         V_S = V[:, S_rel]
#         R_S = Oscar.echelon_form(matrix(F, V_S), reduced = false)
#         any(i -> iszero(R_S[i, i]), 1:(length(S) - 1)) && return false
#         if length(S) <= size(V, 1)
#             !iszero(R_S[length(S), length(S)]) && return false
#         end
#         return true
#     end
#     return false
# end

function is_dual_tropical_root(wd::WalkData, m::MixedCell, d::Vector{Float64})

    d_qq = (QQ).((rationalize.(d)))
    _, w = outer_normal_vector(wd.M, m, d_qq)
    dotps = (i -> dot(w, wd.M.A_Fl[:, i]) + d[i]).(indices(wd.M))
    is_first = true
    last_val = 0.0
    for (j, S) in enumerate(m.inds)
        s1 = first(S)
        for i in m.loc_inds[j]
            i in S && continue
            if dotps[s1] < dotps[i]
                println(m.inds)
                println(i)
                return false
            end
        end
        if !is_first
            if dotps[s1] > last_val
                println(m.inds)
                return false
            end
        end
        last_val = dotps[s1]
        is_first = false
    end
    return true
end

function is_in_mixed_cell_cone(wd::WalkData, m::MixedCell, d::Vector{Float64},
                               p0::DualVector, p1::DualVector)

    tst = true
    for c in keys(wd.walls)
        for (m0, _, sgn) in wd.walls[c]
            if m0 == m
                ineqn_sat = sgn ? dot(d, c) < 0 : dot(d, c) > 0
                if !ineqn_sat
                    println("crossing value t = $(crossing_val(p0, p1, c)), inequality satisfied $(ineqn_sat)")
                end
                tst = tst && ineqn_sat
            end
        end
    end

    return tst 
end
