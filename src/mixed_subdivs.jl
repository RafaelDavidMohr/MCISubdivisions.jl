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
                            ms::Vector{MixedCellInds}, p0::DualVector,
                            p1::DualVector, target_rr_count::Int) where C

    rr_counter = RRCounter(A, V, target_rr_count)

    M = MCI(V, A)
    wd = WalkData(M, [MixedCell(m, M) for m in ms])
    cnt_reached, t_cross = walk_homotopy!(wd, p0, p1, rr_counter = rr_counter)
    return cnt_reached, t_cross, gather_mixed_cells(wd)
end

function starting_system(A::Matrix{Int}, V::Matrix{FqFieldElem})
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

    return homotopy(A, V, A_start, V_start, init_mc, p_start)
end

# works only if conv(A_target) ⊆ conv(A_start)
function homotopy(A_target::Matrix{Int}, V_target::Matrix{FqFieldElem},
                  A_start::Matrix{Int}, V_start::Matrix{FqFieldElem},
                  ms_start::Vector{MixedCellInds},
                  p_start::Vector{Int})

    A_ext = hcat(A_target, A_start)
    V_ext = hcat(V_target, V_start)

    M = MCI(V_ext, A_ext)

    A_size = size(A_target, 2)

    # set up path
    A_ext_size = size(A_ext, 2)
    p0 = forgetful_lift(A_ext_size, collect(1:A_size), p_start)
    p1 = forgetful_lift(A_ext_size, collect(A_size+1:A_ext_size))

    # initial mixed cells
    init_mc = [MixedCell([ci .+ A_size for ci in c], M) for c in ms_start]

    wd = WalkData(M, init_mc)

    return wd, p0, p1
end

# --- Functions related to mixed cell cones --- #

function walk_wall!(wd::WalkData, c::Hyperplane, rr_counter::RRCounter)
    active_mc_data = wd.walls[c]
    delete!(wd.walls, c) 
    new_ms = Set{MixedCellInds}()
    cnt = 0
    for (m, act_index, sgn) in active_mc_data
        delete_mixed_cell!(wd, m, rr_counter)
        mixed_cell_flip!(m, c, wd.M, act_index, sgn, new_ms)
    end
    @info "$(length(new_ms)) new mixed cells"
    for new_mc_inds in new_ms
        new_mc = MixedCell(new_mc_inds, wd.M)
        compute_active_walls!(new_mc, wd.M, wd.walls)
        if rr_counter.target_rr_count >= 0
            add_mixed_cell!(rr_counter, new_mc)
        end
    end
end

function mixed_cell_flip!(m::MixedCell, c::Hyperplane, M::MCI, act_index::Int,
                          sgn::Bool,
                          new_ms::Set{MixedCellInds})

    # indices from which new mixed cell component can come
    # todo: if this doesnt work check if this is correct
    new_inds = setdiff(nz_inds(c), vcat(m.inds...))
    is_exchange = !isempty(new_inds)
    if is_exchange
        @assert isone(length(new_inds))
    else
        new_inds = m.inds[act_index+1]
    end
    
    Ss_new = Vector{Int}[]
    F = prime_field_V(M)
    for (k, i) in enumerate(m.inds[act_index])
        candidate_indices = vcat(m.inds[act_index][1:k-1], m.inds[act_index][k+1:end],
                                 new_inds)
        V_trunc_inds = vcat(candidate_indices,
                            [idx[2:end] for idx in m.inds[1:act_index-1]]...)
        V_trunc = M.V[:, V_trunc_inds]
        K = kernel(matrix(F, V_trunc), side = :right)
        @assert isone(size(K, 2))
        cl = length(candidate_indices)
        S_new = candidate_indices[findall(j -> !iszero(K[j, :]), 1:cl)]
        if partial_sum_sign(S_new, c, sgn)
            push!(Ss_new, S_new)
        end
    end

    inds = vcat(m.inds[act_index], sort(new_inds))
    for S_new in Ss_new
        sort!(S_new)
        S_new_next = sort(setdiff(inds, S_new))
        next_ind = is_exchange ? act_index + 1 : act_index + 2
        new_ms_inds = if length(S_new_next) > 1
            [m.inds[1:act_index-1]..., S_new, S_new_next,
             m.inds[next_ind:end]...]
        else
            [m.inds[1:act_index-1]..., S_new,
             m.inds[next_ind:end]...]
        end
        new_ms_inds in new_ms && continue
        if is_affine_independent(M.A, new_ms_inds)
            push!(new_ms, new_ms_inds)
        end
    end
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

    cayley_config_fac = factorize(cayley_config)
    d = Int(round(det(cayley_config_fac)))
    all_m_inds = vcat(m.inds...)

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
            c = Hyperplane(c_cfs)
            sgn = signbit(partial_sum(m.inds[i], c))
            add_to_dict!(walls, c, (m, i, sgn))
        end
    end
end

# --- Mixed cell data --- #

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

    return outer_normal_vector(M.A, m, d)
end
