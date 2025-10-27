module MCISubdivisions

using Oscar
using AbstractTrees

include("typedefs.jl")
include("helpers.jl")

# --- Main functions --- #

function compute_active_walls_children!(m::MixedCellNode,
                                        circuits_and_walls::Dict{SparseVec, Vector{MixedCellNode}})

    circuits = keys(circuits_and_walls)
    remaining_circuits = filter(c -> iszero(sum(view(c, m.A_indices))), circuits)

    for m_next in m.children
        s_next = m_next.S
        for c in remaining_circuits
            is_active = true
            c_outside_s_next = view(c, setdiff(m.A_indices, s_next))
            signprev = signbit(first(c_outside_s_next))
            for cf in c_outside_s_next[2:end]
                signnext = signbit(cf)
                if signnext != signprev
                    is_active = false
                    break
                end
                signprev = signnext
            end
            is_active && push!(circuits_and_walls[c], m_next)
        end
    end
end

end # module MCISubdivisions
