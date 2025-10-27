struct SparseVec{C}
    cfs::Vector{C}
    inds::Vector{Int}
end

struct WalkData
    A_ext::Matrix{Int}

    mixed_cell_tree::MixedCellNode

    circuit_and_walls::Dict{SparseVec{QQFieldElem}, Vector{MixedCellNode}}
end

function WalkData(A::Matrix{Int},
                  initial_mixed_cell_tree::MixedCellNode)

    A_ext = vcat(A, ones(Int, size(A, 2)))

    @info "computing circuits"
    circuit_indices = circuits(matroid_from_matrix_columns(A_ext))
    circuits = Vector{QQFieldElem}[]
    for c in circuits 
        K = kernel(matrix(QQ, A_ext[:, c], side = :right))
        push!(circuits, SparseVec(K[:, 1], c))
    end
    @info "done, $(length(circuits)) circuits"

    active_walls = Dict([(c, MixedCellNode[]) for c in circuits])
    compute_active_walls_children!(initial_mixed_cell_tree, active_walls)

    return WalkData(A_ext, initial_mixed_cell_tree,
                    circuits, active_walls)
end

struct MixedCellNode 
    S::Vector{Int} # root always with empty S
    A_indices::Vector{Int} # indices into A corresponding to localization at S and its parents
    children::Vector{MixedCellNode}
    active_walls::Vector{Int}
end

AbstractTrees.children(m::MixedCellNode) = m.children
AbstractTrees.childrentype(::Type{<:MixedCellNode}) = Vector{MixedCellNode}
AbstractTrees.ChildIndexing(::Type{<:MixedCellNode}) = AbstractTrees.IndexedChildren()
