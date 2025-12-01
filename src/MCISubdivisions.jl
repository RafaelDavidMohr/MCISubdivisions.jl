module MCISubdivisions

using Oscar
using LinearAlgebra

include("typedefs.jl")
include("helpers.jl")

export mixed_volume, mixed_subdivision

# --- Main functions --- #

function mixed_volume(F::Vector{<:MPolyRingElem})
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    return mixed_volume(A, V)
end

function mixed_subdivision(F::Vector{<:MPolyRingElem})
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    A, mixed_subdivision(A, V)
end

function mixed_volume(A::Matrix{Int}, V::Matrix{C}) where C
    A_ext, cells = mixed_subdivision(A, V)
    return sum([vol(m, A_ext) for m in cells])
end

function mixed_subdivision(A::Matrix{Int}, V::Matrix{C}) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    F = parent(first(Vp))
    rand_mix = matrix(F, (F).(rand(1:characteristic(F)-1, size(V, 1), size(V, 1))))
    Vp = Matrix(rand_mix * matrix(F, Vp))
    A_ext, wd, p0, p1 = total_degree_homotopy(A, Vp)

    walk_homotopy!(wd, p0, p1)

    A_size = size(A, 2)
    result = Set{MixedCell}()
    for c in keys(wd.walls)
        for (m, act_index, sgn) in wd.walls[c]
            any(S -> any(i -> i > A_size, S), m.inds) && continue
            push!(result, m)
        end
    end

    return A_ext, collect(result)
end

function walk_homotopy!(w::WalkData, p0::Vector{Int}, p1::Vector{Int})

    @info "starting homotopy"

    c_prev = nothing
    while true
        c_int = first_intersection(p0, p1, keys(w.walls), c_prev)
        if isnothing(c_int)
            @info "no intersection left, finished"
            return
        end
        @info "intersection found"
        @info "$(length(w.walls[c_int])) mixed cells to flip"
        walk_wall!(w, c_int)
        c_prev = c_int
    end
end

function total_degree_homotopy(A::Matrix{Int}, V::Matrix{FqFieldElem})
    A_size = size(A, 2)

    # extended MCI
    max_deg = maximum(i -> sum(A[:, i]), 1:A_size)
    n = size(A, 1)
    A_ext = copy(A)
    V_ext = copy(V)
    F = parent(first(V))
    for i in 0:n
        A_ext = hcat(A_ext, [j == i ? max_deg : 0 for j in 1:n])
        V_ext = hcat(V_ext, rand_vec_ff(F, n))
    end
    M = MCI(V_ext, A_ext)

    # set up path
    p0 = forgetful_lift(A_size + n + 1, collect(1:A_size))
    p1 = forgetful_lift(A_size + n + 1, collect(A_size+1:A_size+n+1))

    # initial mixed cell
    init_mc = MixedCell([collect(A_size+1:A_size+n+1)], M)
    wd = WalkData(M, [init_mc])

    return A_ext, wd, p0, p1
end

# --- Functions related to mixed cell cones --- #

function walk_wall!(wd::WalkData, c::Hyperplane)
    active_mc_data = wd.walls[c]
    delete!(wd.walls, c) 
    cnt = 0
    for (m, act_index, sgn) in active_mc_data
        delete_mixed_cell!(wd, m)
        new_mixed_cells = mixed_cell_flip(m, c, wd.M, act_index, sgn)
        cnt += length(new_mixed_cells)
        for new_mc in new_mixed_cells
            compute_active_walls!(new_mc, wd.M, wd.walls)
        end
    end
    @info "$(cnt) new mixed cells (possible duplicates)"
end

function mixed_cell_flip(m::MixedCell, c::Hyperplane, M::MCI, act_index::Int, sgn::Bool)

    Mloc = localize(M, m, act_index-1)

    # indices from which new mixed cell component can come
    # todo: if this doesnt work check if this is correct
    new_inds = setdiff(nz_inds(c), vcat(m.inds...))
    is_exchange = !isempty(new_inds)
    inds = if is_exchange 
        sort(union(m.inds[act_index], new_inds))
    else
        sort(vcat(m.inds[act_index], m.inds[act_index+1]))
    end

    Mloc = restrict(Mloc, inds)
    Ss_new = Vector{Int}[]

    for S_new in circuits(Mloc)
        @assert !isone(length(S_new)) "Circuit of length one in localization"
        S_new_rel = Mloc.base_to_rel[S_new]
        sort!(S_new)
        if partial_sum_sign(S_new, c, sgn) && is_affine_independent(Mloc.A_rel, S_new_rel)
            push!(Ss_new, S_new)
        end
    end

    new_mixed_cells = MixedCell[]
    for S_new in Ss_new
        S_new_next = sort(setdiff(inds, S_new))
        next_ind = is_exchange ? act_index + 1 : act_index + 2
        if length(S_new_next) > 1
            new_ms_inds = [m.inds[1:act_index-1]..., S_new, S_new_next,
                           m.inds[next_ind:end]...]
            push!(new_mixed_cells, MixedCell(new_ms_inds, M))
        else
            new_ms_inds = [m.inds[1:act_index-1]..., S_new,
                           m.inds[next_ind:end]...]
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
    F = prime_field_A(M)

    caley_config_modP = Matrix{FqFieldElem}(undef, n + k, 0)
    caley_config_Fl = Matrix{Float64}(undef, n + k, 0)

    shifted_m_indices = Int[]

    projection_to_A = Dict{Int, Int}()

    # build cayley configuration
    end_index = 0
    m_shift_i = Int[]
    for i in 1:k
        apd_modP = [j == i ? one(F) : zero(F) for j in 1:k]
        apd_Fl = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        add_indices = cayley_indices(m, i)
        for (l, j) in enumerate(add_indices)
            if j in m.inds[i]
                push!(m_shift_i, l + end_index)
            end
            projection_to_A[l + end_index] = j
            caley_config_modP = hcat(caley_config_modP,
                                     vcat(M.A_modP[:, j], apd_modP))
            caley_config_Fl = hcat(caley_config_Fl,
                                   vcat(M.A_Fl[:, j], apd_Fl))
        end
        end_index += length(add_indices)
        append!(shifted_m_indices, m_shift_i)
        empty!(m_shift_i)
    end

    # build circuits
    act_index = 1
    shifted_m_index = 1
    for i in 1:size(caley_config_modP, 2)
        if shifted_m_index <= length(shifted_m_indices) && i == shifted_m_indices[shifted_m_index]
            shifted_m_index += 1
            continue
        end
        if iszero(caley_config_modP[n + act_index, i])
            act_index += 1
        end
        circuit_col_indices = vcat(shifted_m_indices[1:shifted_m_index-1], [i],
                                   shifted_m_indices[shifted_m_index:end])
        K_modP = kernel(matrix(F, caley_config_modP[:, circuit_col_indices]),
                        side = :right)
        K_Fl = nullspace(caley_config_Fl[:, circuit_col_indices])
        @assert size(K_modP, 2) == size(K_Fl, 2) == 1 "unexpected dimension in circuit computation"
        K_modP *= (-K_modP[shifted_m_index, 1]^(-1))
        K_Fl *= (-K_Fl[shifted_m_index, 1]^(-1))
        c_cfs_modP = [zero(F) for _ in 1:size(M.A_modP, 2)]
        c_cfs_Fl = zeros(Float64, size(M.A_modP, 2))
        for (j, ind) in enumerate(circuit_col_indices)
            c_cfs_modP[projection_to_A[ind]] += K_modP[j, 1]
            c_cfs_Fl[projection_to_A[ind]] += K_Fl[j, 1]
        end
        c = Hyperplane(c_cfs_modP, c_cfs_Fl)
        sgn = signbit(first(partial_sum(m.inds[act_index], c)))
        add_to_dict!(walls, c, (m, act_index, sgn))
    end
end

# --- MCI functions --- #

function localize(M::RelativeMCI, S::Vector{Int}, rem_inds::Vector{Int})
    S_rel = M.base_to_rel[S]
    rem_inds_rel = M.base_to_rel[rem_inds]
    A = M.A_rel
    V = M.V_rel
    L = linear_span(A, S_rel)
    A_new = project_along_linear_space(L, A[:, rem_inds_rel], length(S_rel) - 1)
    V_new = project_along_linear_space(V[:, S_rel], V[:, rem_inds_rel], length(S_rel) - 1)
    rel_to_base = compose_as_maps(M.rel_to_base, rem_inds_rel)
    return RelativeMCI(M.base, A_new, V_new, rel_to_base)
end

function localize(M::RelativeMCI, m::MixedCell, i::Int)
    M_curr = M
    for j in 1:i
        M_curr = localize(M_curr, m.inds[j], m.loc_inds[j])
    end
    return M_curr
end

function localize(M::MCI, m::MixedCell, i::Int)
    return localize(RelativeMCI(M), m, i)
end

function restrict(M::RelativeMCI, S::Vector{Int})
    S_rel = M.base_to_rel[S]
    return RelativeMCI(M.base, M.A_rel[:, S_rel], M.V_rel[:, S_rel], M.rel_to_base[S_rel])
end

# --- Mixed cell checking/computation for testing --- #

function outer_normal_vector(M::MCI, m::MixedCell,
                             d::Vector{QQFieldElem})

    n = ambient_dim(M)
    F = prime_field_A(M)
    A_modP_lifted = vcat(M.A_modP, transpose((F).(d)))
    A_Fl_lifted = vcat(M.A_Fl, transpose((Float64).(d)))

    eqns_modP = Matrix{FqFieldElem}(undef, n + 1, 0)
    eqns_Fl = Matrix{Float64}(undef, 0, n + 1)
    for S in m.inds
        a_modP = A_modP_lifted[:, first(S)]
        a_Fl = A_Fl_lifted[:, first(S)]
        for i in S[2:end]
            eqns_modP = hcat(eqns_modP, A_modP_lifted[:, i] - a_modP)
            eqns_Fl = vcat(eqns_Fl, transpose(A_Fl_lifted[:, i] - a_Fl))
        end
    end

    K_modP = kernel(matrix(F, eqns_modP))
    @assert isone(size(K_modP, 1)) "mixed cell does not lift to a hyperplane"
    K_modP *= K_modP[1, end]^(-1)

    K_Fl = nullspace(eqns_Fl)
    K_Fl *= K_Fl[end, 1]^(-1)
    return K_modP[1, 1:n], K_Fl[1:n, 1]
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
            
    return MixedCell(result)
end

# checks if m ∪ {S} is a partial mixed cell
function is_partial_mixed_cell(M::MCI, m::MixedCell)

    M_curr = RelativeMCI(M)

    for j in 1:length(m.inds)
        !is_partial_mixed_cell(M_curr, m.inds[j]) && return false
        M_curr = localize(M_curr, m.inds[j], m.outside_inds[j])
    end

    return true
end

function is_partial_mixed_cell(M::RelativeMCI,
                               S::Vector{Int})

    S_rel = M.base_to_rel[S]
    A = M.A_rel
    V = M.V_rel
    if is_affine_independent(A, S_rel)
        F = parent(first(V))
        V_S = V[:, S_rel]
        R_S = Oscar.echelon_form(matrix(F, V_S), reduced = false)
        any(i -> iszero(R_S[i, i]), 1:(length(S) - 1)) && return false
        if length(S) <= size(V, 1)
            !iszero(R_S[length(S), length(S)]) && return false
        end
        return true
    end
    return false
end

function is_in_mixed_cell_cone(wd::WalkData, m::MixedCell, d::Vector{Float64},
                               p0::Vector{Int}, p1::Vector{Int})
    d_qq = (QQ).((rationalize).(d))
    _, w = outer_normal_vector(wd.M, m, d_qq)
    sinds = sort(indices(wd.M), by = i -> dot(w, wd.M.A_Fl[:, i]) + d[i], rev = true)
    ml = sum((length).(m.inds))
    # println(sinds[1:ml])
    # println(m)
    # println(m)
    # for i in indices(M)
    tst = true
    k = 0
    for S in m.inds
        tst = tst && S == sort(sinds[1+k:length(S)+k])
        k += length(S)
    end

    if !tst
    # println("----")
        tst2 = true
        for c in keys(wd.walls)
            for (m0, _, sgn) in wd.walls[c]
                if m0 == m
                    println("crossing value t = $(crossing_val(p0, p1, c))")
                    tst2 = tst2 && (sgn ? dot(d, c) > 0 : dot(d, c) < 0)
                end
            end
        end
        if tst2
            @info "inconsistent test result"
        end
        return false
    end
    return true
end

end # module MCISubdivisions
