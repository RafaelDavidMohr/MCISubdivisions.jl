module MCISubdivisions

using Oscar
using LinearAlgebra

include("typedefs.jl")
include("helpers.jl")

# --- Main functions --- #

# --- Functions related to mixed cell cones --- #

function compute_active_walls!(m::MixedCell,
                               M::MCI,
                               walls::Dict{Circuit, Vector{Tuple{Int, MixedCell}}})

    k = length(m)
    n = ambient_dim(M)
    F = base_ring(first(M.A_modP))

    caley_config_modP = Matrix{FqFieldElem}(undef, 0, n + k)
    caley_config_Fl = Matrix{Float64}(undef, 0, n + k)

    remaining_indices = collect(1:size(M.A_modP, 2))
    shifted_m_indices = Vector{Int}[]

    projection_to_A = Dict{Int, Int}()

    # build cayley configuration
    end_index = 0
    m_shift_i = Int[]
    for i in 1:k
        apd_modP = [j == i ? one(F) : zero(F) for j in 1:k]
        apd_Fl = [j == i ? one(Float64) : zero(Float64) for j in 1:k]
        for (l, j) in enumerate(remaining_indices)
            if j in m.inds[i]
                push!(m_shift_i, l + end_index)
            end
            projection_to_A[j + end_index] = l
            caley_config_modP = hcat(caley_config_modP, vcat(M.A_modP[:, j], apd_modP))
            caley_config_Fl = hcat(caley_config_Fl, vcat(M.A_Fl[:, j], apd_Fl))
        end
        end_index += length(remaining_indices)
        setdiff!(remaining_indices, m.inds[i]) # CAREFUL
        push!(shifted_m_indices, copy(m_shift_i))
    end

    # build circuits
    m_index = 1
    shifted_m_index = 1
    for i in 1:size(caley_config_modP, 2)
        if i == shifted_m_indices[shifted_m_index]
            shifted_m_index += 1
            continue
        end
        if iszero(caley_config_modP[i, n + m_index])
            m_index += 1
        end
        circuit_col_indices = vcat(shifted_m_indices[1:shifted_m_index-1], [i], shifted_m_indices[shifted_m_index:end])
        K_modP = kernel(matrix(F, caley_config_modP[:, circuit_col_indices]), side = :right)
        K_Fl = nullspace(caley_config_Fl[:, circuit_col_indices])
        @assert size(K_modP, 2) == size(K_Fl, 2) == 1 "unexpected dimension in circuit computation"
        c_cfs_modP = [zero(F) for _ in 1:size(M.A_modP, 2)]
        c_cfs_Fl = zeros(Float64, size(M.A_modP, 2))
        for (j, ind) in enumerate(circuit_col_indices)
            c_cfs_modP[projection_to_A[ind]] += K_modP[j, 1]
            c_cfs_Fl[projection_to_A[ind]] += K_Fl[j, 1]
        end
        c_final_inds = findall(!iszero, c_cfs_modP)
        c = Circuit(c_final_inds, c_cfs_modP[c_final_inds], c_cfs_Fl[c_final_inds])
        add_to_dict!(walls, c, (m_index, m))
    end
end

# --- Mixed cell checking/computation --- #

# optimization: output one floating point, one finite field representation
function outer_normal_vector(M::MCI, m::MixedCell,
                             d::Vector{QQFieldElem})

    n = ambient_dim(M)
    F = base_ring(first(M.A_modP))
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

# this should work for general d, assuming that w is a tropical root?
function find_dual_tropical_root(M::MCI, d::Vector{QQFieldElem},
                                 w::Vector{QQFieldElem})

    result = Vector{Int}[]

    A_card = size(M.A_ext, 2)
    n = size(M.A_ext, 1) - 1
    s = sort(1:A_card, by = i -> dot(w, M.A_ext[1:n, i]) + d[i], rev = true)
    println(s)

    deg = dot(w, M.A_ext[1:n, first(s)]) + d[first(s)]
    codim = 0
    Sj = Int[]

    i = 1
    while codim < n
        if i <= A_card && dot(w, M.A_ext[1:n, s[i]]) + d[s[i]] == deg
            push!(Sj, s[i])
        else
            push!(result, copy(Sj))
            codim += length(Sj) - 1
            Sj = Int[]
            if i < A_card
                deg = dot(w, M.A_ext[1:n, s[i+1]]) + d[s[i+1]]
            end
        end
        i += 1
    end
        
    return result
end

function is_partial_mixed_cell(M::MCI, prnt::MixedCellNode, S::Vector{Int})

    ancestor_cell_indices = Vector{Int}[]
    node = prnt
    while !isnothing(node)
        !isempty(node.S) && pushfirst!(ancestor_cell_indices, node.S)
        node = node.parent
    end

    # check affine independence
    F, A_extF = reduce_mod_rand_prime(M.A_ext)
    n = size(M.V, 1)
    expected_dimension = sum((length).(ancestor_cell_indices)) - length(ancestor_cell_indices)
    expected_dimension += length(S) - 1
    Oscar.rank(matrix(F, M.A_ext[:, vcat(ancestor_cell_indices..., S)])) - 1 != expected_dimension && return false

    # check rank condition on matroid
    F = parent(first(M.V))
    V_S = M.V[:, vcat(S, ancestor_cell_indices...)]
    R_S = Oscar.echelon_form(matrix(F, V_S))
    any(i -> iszero(R_S[i, i]) || iszero(R_S[i, length(S)]), 1:(length(S) - 1)) && return false

    return true
end

end # module MCISubdivisions
