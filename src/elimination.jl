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

function symbolic_volume(A_mod::Matrix{C},
                         m::Vector{Vector{Int}},
                         evl::Vector{Int}) where C

    mat = hcat([linear_span(A_mod, S) for S in m]...)
    R = parent(first(A_mod))
    d = det(matrix(R, mat))
    return d(evl...) > 0 ? d : -d
end

function elim_supp_func(E::ElimData)
    R, t = polynomial_ring(QQ, "t")
    n = size(E.V, 1) - 1
    Aw = vcat((R).(E.A[1:n, :]),
              t .* permutedims(vcat(zeros(Int, n), E.current_covec)) * E.A)
    AC = vcat((R).(E.A[1:n, :]), (R).(repeat([E.shift], 1, size(E.A, 2))))
    A_mod = hcat(Aw, AC)
    sv = sum([symbolic_volume(A_mod, m, [1]) for m in E.current_ms])
    return Int(coeff(sv, 1))
end

function elim_vertex(E::ElimData)
    n = size(E.V, 1) - 1
    k = size(E.A, 1) - n - 1
    R, t = polynomial_ring(QQ, ["t$i" for i in 1:k+1])
    Aw = vcat((R).(E.A[1:n, :]),
              permutedims(vcat(zeros(R, n), t)) * (R).(E.A))
    AC = vcat((R).(E.A[1:n, :]), (R).(repeat([E.shift], 1, size(E.A, 2))))
    A_mod = hcat(Aw, AC)
    sv = sum([symbolic_volume(A_mod, m, E.current_covec) for m in E.current_ms])
    return [Int(coeff(sv, v)) for v in t]
end
    
