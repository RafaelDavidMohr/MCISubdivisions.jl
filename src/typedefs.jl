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

mutable struct MixedCellNode 
    S::Vector{Int} # root always with empty S

    A_remaining::Vector{Int} # indices corresponding to support of localization at S and its parents
    circuits::Vector{SparseVec{QQFieldElem}}
    circuit_indices::Vector{Int} # indices of circuits which give circuits of localization at S and its parents
    
    parent::Union{Nothing, MixedCellNode}
    children::Vector{MixedCellNode}
end

AbstractTrees.children(m::MixedCellNode) = m.children
AbstractTrees.childrentype(::Type{<:MixedCellNode}) = Vector{MixedCellNode}
AbstractTrees.ChildIndexing(::Type{<:MixedCellNode}) = AbstractTrees.IndexedChildren()
AbstractTrees.ParentLinks(::Type{<:MixedCellNode}) = AbstractTrees.StoredParents()

struct WalkData
    M::MCI

    mixed_cell_tree::MixedCellNode

    walls::Dict{Vector{QQFieldElem}, Vector{MixedCellNode}}
end

function WalkData(M::MCI,
                  initial_mixed_cell::Vector{Vector{Int}})

    A_ext = M.A_ext

    # compute affine circuits of A
    @info "computing affine circuits"
    F, A_extF = reduce_mod_rand_prime(A_ext)
    circuit_indices = circuits(matroid_from_matrix_columns(matrix(F, A_extF)))
    circs = Vector{QQFieldElem}[]
    for c in circs 
        K = kernel(matrix(QQ, A_ext[:, c], side = :right))
        push!(circs, SparseVec(K[:, 1], c))
    end
    @info "done, $(length(circs)) affine circuits"

    # set up tree encoding initial mixed cell
    root = MixedCellNode(Int[], collect(1:size(A_ext, 2)), circs,
                         collect(1:length(circs)),
                         nothing, MixedCellNode[])
    node = root
    for S in initial_mixed_cell
        new_circuit_indices = filter(i -> iszero(sum(subvec(circs[i], S))),
                                     node.circuit_indices)
        new_node = MixedCellNode(S, setdiff(node.A_remaining, S),
                                 circs, new_circuit_indices,
                                 node, MixedCellNode[])
        node.children = [new_node]
        node = new_node
    end

    # compute set of active wall at initial mixed cell
    walls = Dict{Vector{QQFieldElem}, Vector{MixedCellNode}}()
    for m in AbstractTrees.PreOrderDFS(root)
        compute_active_walls_children!(m, walls, size(A_ext, 2))
    end

    return WalkData(M, root, walls)
end
