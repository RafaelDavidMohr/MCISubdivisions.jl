function mixed_subdivision_direct(A::Matrix{Int}, V::Matrix{C}, d::Vector{Int}) where C
    Ad = vcat(A, permutedims(d))
    @info "computing upper facets"
    P = convex_hull(transpose(Ad))
    fs = collect(facets(Halfspace, P))
    filter!(f -> f.a[1, end] > 0, fs)
    fs_points = Vector{Int}[]
    for f in fs
        dts = vec(Matrix(f.a) * Ad)
        mv = maximum(dts)
        fs_inds = findall(i -> dts[i] == mv, 1:size(A, 2))
        if length(fs_inds) == size(A, 1) + 1 
            push!(fs_points, fs_inds)
        end
    end

    F = parent(first(V))
    init_comps = unique([f[findall(!iszero, kernel(matrix(F, V[:, f]), side = :right))] for f in fs_points])
    filter!(f -> length(f) != size(A, 1), init_comps)
    return init_comps
end
