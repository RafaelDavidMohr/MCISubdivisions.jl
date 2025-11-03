module MCISubdivisions

using Oscar
using AbstractTrees

include("typedefs.jl")
include("helpers.jl")

# --- Main functions --- #

# --- Functions related to mixed cell cones --- #

function compute_active_walls_children!(m::MixedCellNode,
                                        walls::Dict{Vector{QQFieldElem}, Vector{MixedCellNode}},
                                        A_length::Int)

    for m_next in m.children
        s_next = m_next.S
        for c in view(m.circuits, m.circuit_indices)
            is_active = true
            c_outside_s_next = subvec(c, setdiff(m.A_remaining, s_next))
            signprev = signbit(first(c_outside_s_next))
            for cf in c_outside_s_next[2:end]
                signnext = signbit(cf)
                if signnext != signprev
                    is_active = false
                    break
                end
                signprev = signnext
            end
            if is_active
                wall_vector = densify(c[m.A_remaining], A_length)
                add_to_dict!(walls, signprev ? wall_vector : -wall_vector, m_next)
            end
        end
    end
end

# --- Mixed cell checking/computation --- #

# this should work for general d, assuming that w is a tropical root?
function find_dual_tropical_root(M::MCI, d::Vector{QQFieldElem},
                                 w::Vector{QQFieldElem})

    result = Vector{Int}[]

    A_card = size(M.A_ext, 2)
    n = size(M.A_ext, 1) - 1
    s = sort(1:A_card, by = i -> dot(w, M.A_ext[1:n, i]) + d[i], rev = true)

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

function is_partial_mixed_cell(M::MCI, parent::MixedCellNode, S::Vector{Int})

    ancestor_cell_indices = Vector{Int}[]
    node = parent
    while !isnothing(node)
        pushfirst!(ancestor_cell_indices, node.S)
        node = node.parent
    end

    # check affine independence
    n = size(M.V, 1)
    expected_dimension = sum((length).(ancestor_cell_indices)) - length(ancestor_cell_indices)
    expected_dimension += length(S) - 1
    rank(M.A_ext[:, vcat(ancestor_cell_indices..., S)]) - 1 != expected_dimension && return false

    # check rank condition on matroid
    V_S = M.V[:, vcat(S, ancestor_cell_indices...)]
    R_S = reduced_echelon_form(V_S)
    any(i -> iszero(R_S[i, i]) || iszero(R_S[i, length(S)]), 1:length(S)) && return false

    return true
end

end # module MCISubdivisions
