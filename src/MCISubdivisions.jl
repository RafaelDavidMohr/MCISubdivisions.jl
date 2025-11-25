module MCISubdivisions

using Oscar
using LinearAlgebra

include("typedefs.jl")
include("helpers.jl")

# --- Main functions --- #

function mixed_subdivision(A::Matrix{Int}, V::Matrix{C}) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    F = parent(first(V))
    rand_mix = matrix(F, (F).(rand(1:characteristic(F)-1, size(V, 1), size(V, 1))))
    Vp = rand_mix * Vp
    wd, path = total_degree_homotopy(A, Vp)
    walk_homotopy!(wd, path)
    return unique!(vcat([(tup -> tup[1]).(wd.walls[c]) for c in keys(wd.walls)]...))
end

function walk_homotopy!(w::WalkData, p::HomotopyPath)

    @info "starting homotopy"

    while !is_completed(path)
        c_int = first_intersection_with_path!(p, keys(w.walls))
        if isnothing(c_int)
            @info "no intersection left, finished"
            return
        end
        @info "intersection found, $(length(p.points) - 1) segments remaining"
        @info "$(length(w.walls[c_int])) mixed cells to flip"
        walk_wall!(w, c_int)
    end
end

function total_degree_homotopy(A::Matrix{Int}, V::Matrix{FqFieldElem})
    A_size = size(A, 2)
    zeropos = findfirst(i -> iszero(A[:, i]), 1:A_size)
    @assert !isnothing(zeropos) "Support does not contain origin, consider shifting"

    # extended MCI
    max_deg = max(i -> sum(A[:, i]), 1:A_size)
    n = size(A, 1)
    A_ext = copy(A)
    V_ext = copy(V)
    F = parent(first(V))
    for i in 1:n
        A_ext = hcat(A_ext, [j == i ? 1 : 0 for j in 1:n])
        V_ext = hcat(V_ext, (F).(rand(1:characteristic(F)-1, n)))
    end
    M = MCI(V_ext, A_ext)

    # set up path
    p0 = vcat(zeros(Float64, zeropos - 1),
              [100 * rand()],
              zeros(Float64, A_size - zeropos),
              100 * rand(n))
    p1 = vcat(100 * rand(A_size), zeros(Float64, n))
    path = piecewice_linear_path(p0, p1)

    # initial mixed cell
    init_mc = MixedCell([zeropos, collect(A_size+1:A_size+1+n)...])
    wd = WalkData(M, [init_mc])

    return wd, path
end

# --- Functions related to mixed cell cones --- #

function walk_wall!(wd::WalkData, c::Circuit)
    active_mc_data = wd.walls[c]
    delete!(wd.walls, c)
    cnt = 0
    for (m, act_index, sgn) in active_mc_data
        new_mixed_cells = mixed_cell_flip(m, c, act_index, sgn)
        cnt += length(new_mixed_cells)
        for new_mc in new_mixed_cells
            compute_active_walls!(new_mc, wd.M, wd.walls)
        end
    end
    @info "$(cnt) new mixed cells (possible duplicates)"
end

function mixed_cell_flip(m::MixedCell, c::Circuit, M::MCI, act_index::Int, sgn::Bool)

    A_loc, V_loc, rem_inds, ind_map = localize(M.A_modP, M.V, m.inds[1:act_index-1])

    inds = if act_index == length(m)
        excld = isone(act_index) ? Int[] : vcat(m.inds[1:act_index-1]...)
        setdiff(union(m.inds[act_index], c.inds), excld)
    else
        vcat(m.inds[act_index], m.inds[act_index + 1])
    end

    matr = matroid_from_matrix_columns(matrix(prime_field_V(M), V_loc))
    matr = restriction(matr, [ind_map[i] for i in inds])

    Ss_new = Vector{Int}[]

    for S_new in circuits(matr)
        S_new_A_inds = rem_inds[S_new]
        if partial_sum_sign(S_new_A_inds, c, sgn) && is_affine_independent(A_loc, S_new)
            push!(Ss_new, S_new_A_inds)
        end
    end

    new_mixed_cells = MixedCell[]
    for S_new in Ss_new
        S_new_next = setdiff(inds, S_new)
        if length(S_new_next) > 1
            push!(new_mixed_cells, MixedCell([m.inds[1:act_index-1]..., S_new, S_new_next, m.inds[act_index + 2:end]...]))
        else
            push!(new_mixed_cells, MixedCell([m.inds[1:act_index-1]..., S_new, m.inds[act_index + 2:end]...]))
        end
    end

    return new_mixed_cells
end

function compute_active_walls!(m::MixedCell,
                               M::MCI,
                               walls::Dict{Circuit, Set{Tuple{MixedCell, Int, Bool}}})

    k = length(m)
    n = ambient_dim(M)
    F = prime_field_A(M)

    caley_config_modP = Matrix{FqFieldElem}(undef, n + k, 0)
    caley_config_Fl = Matrix{Float64}(undef, n + k, 0)

    shifted_m_indices = Int[]

    projection_to_A = Dict{Int, Int}()

    # build cayley configuration
    # Here: to the ith part of the Cayley configuration we just need to add
    # the first element of S_i+1. To the last part of the Cayley configuration
    # we add all the elements that lie outside of S1 u ... u Sk
    # TODO: this loop can probably be optimized
    end_index = 0
    m_shift_i = Int[]
    for i in 1:k
        apd_modP = [j == i ? one(F) : zero(F) for j in 1:k]
        apd_Fl = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        add_indices = if i < k
            vcat(m.inds[i], [first(m.inds[i+1])])
        else
            vcat(m.inds[i], setdiff(collect(1:size(M.A_modP, 2)), vcat(m.inds...)))
        end
        for (l, j) in enumerate(add_indices)
            if j in m.inds[i]
                push!(m_shift_i, l + end_index)
            end
            projection_to_A[l + end_index] = j
            caley_config_modP = hcat(caley_config_modP, vcat(M.A_modP[:, j], apd_modP))
            caley_config_Fl = hcat(caley_config_Fl, vcat(M.A_Fl[:, j], apd_Fl))
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
        circuit_col_indices = vcat(shifted_m_indices[1:shifted_m_index-1], [i], shifted_m_indices[shifted_m_index:end])
        K_modP = kernel(matrix(F, caley_config_modP[:, circuit_col_indices]), side = :right)
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
        c = Circuit(c_cfs_modP, c_cfs_Fl)
        sgn = signbit(first(partial_sum(m.inds[act_index], c)))
        add_to_dict!(walls, c, (m, act_index, sgn))
    end
end

# --- MCI functions --- #

function localize(A::Matrix{FqFieldElem}, V::Matrix{FqFieldElem}, S::Vector{Int})
    L = linear_span(A, S)
    rem_inds = setdiff(collect(1:size(A, 2)), S)
    A_new = project_along_linear_space(L, A[:, rem_inds], length(S) - 1)
    V_new = project_along_linear_space(V[:, S], V[:, rem_inds], length(S) - 1)
    ind_map = Dict{Int, Int}()
    for (i, j) in enumerate(rem_inds)
        ind_map[j] = i
    end
    return A_new, V_new, rem_inds, ind_map
end

function localize(A::Matrix{FqFieldElem}, V::Matrix{FqFieldElem}, m::Vector{Vector{Int}})
    A_curr, V_curr, rem_inds = A, V, collect(1:size(A,2))
    for S in m
        A_curr, V_curr, rem_inds, _ = localize(A_curr, V_curr, S)
    end
    ind_map = Dict{Int, Int}()
    for (i, j) in enumerate(rem_inds)
        ind_map[j] = i
    end
    return A_curr, V_curr, rem_inds, ind_map
end

# --- Mixed cell checking/computation for testing --- #

function outer_normal_vector(M::MCI, m::MixedCell,
                             d::Vector{QQFieldElem})

    n = ambient_dim(M)
    F = prime_field_A(M)
    A_modP_lifted = vcat(M.A_modP, (F).(d))
    A_Fl_lifted = vcat(M.A_Fl, (Float64).(d))

    eqns_modP = Matrix{FqFieldElem}(undef, n + 1, 0)
    eqns_Fl = Matrix{Float64}(undef, 0, n + 1)
    for S in m.inds
        a_modP = A_modP_lifted[:, first(S)]
        a_Fl = A_Fl_lifted[:, first(S)]
        for i in S[2:end]
            eqns_modP = hcat(eqns, A_modP_lifted[:, i] - a_modP)
            eqns_Fl = vcat(eqns, transpose(A_Fl_lifted[:, i] - a_Fl))
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

    inds = copy(m.inds)
    A, V = M.A_modP, M.V

    for j in 1:length(inds)
        !is_partial_mixed_cell(A, V, inds[j]) && return false
        A, V, _, ind_map = localize(A, V, inds[j])
        for l in j+1:length(inds)
            inds[l] = [ind_map[a] for a in inds[l]]
        end
    end

    return true
end

function is_partial_mixed_cell(A::Matrix{FqFieldElem},
                               V::Matrix{FqFieldElem},
                               S::Vector{Int})

    if is_affine_independent(A, S)
        F = parent(first(V))
        V_S = V[:, S]
        R_S = Oscar.echelon_form(matrix(F, V_S), reduced = false)
        any(i -> iszero(R_S[i, i]), 1:(length(S) - 1)) && return false
        if length(S) <= size(V, 1)
            !iszero(R_S[length(S), length(S)]) && return false
        end
        return true
    end
    return false
end

end # module MCISubdivisions
