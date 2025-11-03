struct SparseVec{C}
    cfs::Vector{C}
    inds::Vector{Int}
end

struct MCI
    V::Matrix{FqFieldElem} # stored over random finite field to speed up computations
    A_ext::Matrix{Int64}

    function MCI(V::Matrix{FqFieldElem}, A::Matrix{Int64})
        A_ext = vcat(A, ones(Int64, 1, size(A, 2)))
        return new(V, A_ext)
    end
end

function MCI(V::Matrix{QQFieldElem}, A::Matrix{Int64})
    _, Vp = reduce_mod_rand_prime(V)
    return MCI(Vp, A)
end

function WalkData(A::Matrix{Int},
                  initial_mixed_cell::Vector{Vector{Int}})

    A_ext = vcat(A, ones(Int, size(A, 2)))

    # compute affine circuits of A
    @info "computing affine circuits"
    circuit_indices = circuits(matroid_from_matrix_columns(A_ext))
    circuits = Vector{QQFieldElem}[]
    for c in circuits 
        K = kernel(matrix(QQ, A_ext[:, c], side = :right))
        push!(circuits, SparseVec(K[:, 1], c))
    end
    @info "done, $(length(circuits)) affine circuits"

    # set up tree encoding initial mixed cell
    root = MixedCellNode(Int[], collect(1:size(A, 2)), circuits,
                         collect(1:length(circuits)),
                         nothing, MixedCellNode[])
    node = root
    for S in initial_mixed_cell
        new_circuit_indices = filter(i -> iszero(sum(subvec(circuits[i], S))),
                                     node.circuit_indices)
        new_node = MixedCellNode(S, setdiff(node.A_remaining, S),
                                 circuits, new_circuit_indices,
                                 node, MixedCellNode[])
        node.children = [new_node]
        node = new_node
    end

    # compute set of active wall at initial mixed cell
    walls = Dict{Vector{QQFieldElem}, Vector{MixedCellNode}}()
    for m in AbstractTrees.PreOrderDFS(root)
        compute_active_walls_children!(m, walls)
    end

    return WalkData(A_ext, root, walls)
end

mutable struct MixedCellNode 
    S::Vector{Int} # root always with empty S

    A_remaining::Vector{Int} # indices corresponding to support of localization at S and its parents
    circuits::Vector{SparseVec{QQFieldElem}}
    circuit_indices::Vector{Int} # indices of circuits which give circuits of localization at S and its parents
    
    parent::Union{Nothing, MixedCellNode}
    children::Vector{MixedCellNode}
end

struct WalkData
    A_ext::Matrix{Int64}

    mixed_cell_tree::MixedCellNode

    walls::Dict{Vector{QQFieldElem}, Vector{MixedCellNode}}
end

AbstractTrees.children(m::MixedCellNode) = m.children
AbstractTrees.childrentype(::Type{<:MixedCellNode}) = Vector{MixedCellNode}
AbstractTrees.ChildIndexing(::Type{<:MixedCellNode}) = AbstractTrees.IndexedChildren()
AbstractTrees.ParentLinks(::Type{<:MixedCellNode}) = AbstractTrees.StoredParents()
