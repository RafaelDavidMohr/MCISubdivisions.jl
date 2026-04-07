# --- Functions for tropical elimination --- #

# function construct_polytope!(E::ElimData)

#     n = size(E.M.V, 1) - 1 # n + 1 input equations
#     k = size(E.A, 1) - n - 1
#     amb_dim = k + 1 
#     @info "computing initial vertices"

#     w = rand(-1000:1000, amb_dim)
#     P = convex_hull([elim_vertex!(E, w)])
#     dm = 0
#     while dm < amb_dim 
#         @info "dimension $(dm)"
#         af = affine_hull(P)
#         cfs = rand(-1000:1000, length(af))
#         w = sum(cfs .* [intify(h.a[1, :]) for h in af])
#         vert = elim_vertex!(E, w)
#         if all(h -> vert in h, af) 
#             vert = elim_vertex!(E, -w)
#             all(h -> vert in h, af) && break
#         end
#         P = convex_hull(P, convex_hull([vert]))
#         dm = dim(P)
#     end
#     @info "done"

#     facts_confirmed = AffineHalfspace{QQFieldElem}[]

#     all_confirmed = false
#     while !all_confirmed
#         facts = facets(P)
#         all_confirmed = true

#         for fc in facts
#             fc in facts_confirmed && continue

#             nv = (Int).(fc.a[1,:])
#             val = fc.b
#             val2 = elim_supp_func!(E, nv)
#             if val2 == val
#                 @info "facet confirmed"
#                 push!(facts_confirmed, fc)
#                 continue
#             end

#             new_vert = elim_vertex!(E, nv)

#             if !(new_vert in P) # check if new vertex was actually obtained
#                 @info "new vertex"
#                 P = convex_hull(P, convex_hull([new_vert]))
#                 all_confirmed = false
#                 break
#             end
#         end
#     end

#     return P
# end

function elim_supp_func(A::Matrix{Int}, V::Matrix{C}, w::Vector{Int}) where C
    n = size(V, 1) - 1
    Aw = vcat(A[1:n, :], permutedims(vcat(zeros(Int, n), w)) * A)
    return mixed_shadow_volume(Aw, V)
end

# mixed shadow volume w.r.t. last coordinate
function mixed_shadow_volume(A::Matrix{Int}, V::Matrix{C}, cashed_mv::Int, cashed_shft::Int) where C
    A_proj = copy(A)
    A_proj[end, :] = zeros(Int, 1, size(A, 2))
    min_exp = minimum(i -> A[end, i], 1:size(A, 2))

    with_logger(NullLogger()) do
        if min_exp >= 0
            return mixed_volume(hcat(A, A_proj), hcat(V, V))
        end

        if min_exp < cashed_shft
            error("compute new shift")
        end

        A_proj_shft = copy(A)
        A_proj_shft[end, :] = repeat([cashed_shft], 1, size(A, 2))

        a = mixed_volume(hcat(A, A_proj_shft), hcat(V, V))
        # b = mixed_volume(hcat(A_proj, A_proj_shft), hcat(V, V))
        return a - cashed_mv
    end
end
