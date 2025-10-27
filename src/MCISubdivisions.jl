module MCISubdivisions

using Oscar
using AbstractTrees

include("typedefs.jl")
include("helpers.jl")

# --- Main functions --- #

function compute_active_walls_children!(M::MixedCellNode,
                                        circuits_and_walls::Dict{SparseVec, Vector{MixedCellNode}})

    remaining_circuits = filter(c -> iszero(sum(view(c, M.A_indices))), circuits)

    for M_next in M.children
        S_next = M_next.S
        for c in remaining_circuits
            is_active = true
            c_outside_S_next = view(c, setdiff(M.A_indices, S_next))
            signprev = signbit(first(c_outside_S_next))
            for cf in c_outside_S_next[2:end]
                signnext = signbit(cf)
                if signnext != signprev
                    is_active = false
                    break
                end
                signprev = signnext
            end
            is_active && push!(circuits_and_walls[c], M_next)
        end
    end
end

end # module MCISubdivisions
