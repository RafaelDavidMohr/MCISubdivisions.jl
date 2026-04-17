function mixed_subdivision_direct(A::Matrix{Int}, V::Matrix{C}, d::Vector{Int}) where C
    @info "dimension $(size(A, 1)), size $(size(A, 2))"
    comps = mixed_cell_components(A, V, d)
    @info "$(length(comps)) candidates"
    rslt = Vector{Vector{Int}}[]
    for c in comps
        if length(c) == size(A, 1) + 1
            push!(rslt, [c])
            continue
        end
        rem_inds = setdiff(1:size(A, 2), c)
        A_loc, V_loc = localize(A, V, c, rem_inds)
        d_loc = d[rem_inds]
        ms_loc = mixed_subdivision_direct(A_loc, V_loc, d_loc)
        isempty(ms_loc) && continue
        for c_rem in ms_loc
            full_cell = vcat([c], [rem_inds[cr] for cr in c_rem])
            full_cell_inds = sort(vcat(full_cell...))
            fcl = length(full_cell_inds)
            w = outer_normal_vector(A, full_cell, d)
            sp = sortperm([dot(A[:, i], w) + d[i] for i in 1:size(A, 2)], rev = true)
            if full_cell_inds == sort(sp[1:fcl])
                push!(rslt, full_cell)
            end
        end
    end
    @info "$(length(rslt)) mixed cells"
    return rslt
end

function mixed_cell_components(A::Matrix{Int}, V::Matrix{C}, d::Vector{Int}) where C
    col_inds = select_max_weight_columns(A, d)
    sd = subdivision_of_points(transpose(A[:, col_inds]), -d[col_inds])
    cls = [col_inds[c] for c in maximal_cells(sd)]

    F = parent(first(V))
    init_comps = Vector{Int}[]
    for c in cls
        k = kernel(matrix(F, V[:, c]), side = :right)[:, 1]
        push!(init_comps, sort(c[findall(!iszero, k)]))
    end
    filter!(f -> length(f) != 1, init_comps)
    return unique(init_comps)
end

function mixed_cell_components(A::Matrix{Int}, V::Matrix{C}, d::Vector{Int}, cell_comp::Vector{Int}) where C
    rem_inds = setdiff(1:size(A, 2), cell_comp)
    A_loc, V_loc = localize(A, V, cell_comp, rem_inds)
    d_loc = d[rem_inds]
    return [rem_inds[c] for c in mixed_cell_components(A_loc, V_loc, d_loc)]
end

function localize(A::Matrix{Int}, V::Matrix{C}, cell_comp::Vector{Int}, rem_inds::Vector{Int}) where C

    
    a0 = A[:, first(cell_comp)]
    A_shft = hcat([A[:, i] - a0 for i in 1:size(A, 2)]...)
    A_shft = A_shft[:, vcat(cell_comp, rem_inds)]
    A_shft_hf = Matrix{Int}(hermite_form(matrix(ZZ, A_shft)))

    dm = length(cell_comp) - 1
    A_loc = A_shft_hf[dm+1:end, dm+2:end]

    F = parent(first(V))
    V_perm = V[:, vcat(cell_comp, rem_inds)]
    V_perm_re = Matrix(echelon_form(matrix(F, V_perm)))
    display(V_perm_re)
    V_loc = V_perm_re[dm+1:end, dm+2:end]

    return A_loc, V_loc
end
