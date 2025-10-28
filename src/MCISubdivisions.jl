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

# this should work for general d?
function find_dual_tropical_root(M::MCI, d::Vector{QQFieldElem},
                                 w::Vector{QQFieldElem})

    result = Vector{Int}[]

    A_card = size(M.A_ext, 2)
    n = size(M.A_ext, 1) - 1
    s = sort(1:A_card, by = i -> dot(w, M.A_ext[1:n, i]) + d[i], rev = true)

    d = dot(w, M.A_ext[1:n, first(s)]) + d[first(s)]
    codim = 0
    Sj = Int[]

    i = 1
    while codim < n
        if dot(w, M.A_ext[1:n, s[i]]) + d[s[i]] == d
            push!(Sj, s[i])
        else
            push!(result, copy(Sj))
            codim += length(Sj) - 1
            Sj = Int[]
        end
        i += 1
    end
        
    return result
end
    
    
                                 

end # module MCISubdivisions
