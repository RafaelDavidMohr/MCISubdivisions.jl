# --- Functions for tropical elimination --- #

function construct_polytope!(E::ElimData)

    n = size(E.V, 1) - 1 # n + 1 input equations
    k = size(E.A, 1) - n - 1
    amb_dim = k + 1 
    @info "computing initial vertices"

    w = rand(-10:10, amb_dim)
    P = convex_hull([elim_vertex(E)])
    dm = 0
    while dm < amb_dim 
        @info "dimension $(dm)"
        af = affine_hull_int(P)
        cfs = rand(-10:10, length(af))
        w = (Int).(sum(cfs .* [(numerator).(h.a[1, :]) for h in af]))
        vert = elim_vertex!(E, w)
        if all(h -> vert in h, af)
            vert = elim_vertex!(E, -w)
            all(h -> vert in h, af) && break
        end
        P = convex_hull(P, convex_hull([vert]))
        dm = dim(P)
    end
    @info "done"

    facts_confirmed = AffineHalfspace{QQFieldElem}[]

    all_confirmed = false
    while !all_confirmed
        facts = facets(P)
        all_confirmed = true

        for fc in facts
            fc in facts_confirmed && continue

            nv = (Int).(fc.a[1,:])
            val = fc.b
            val2 = elim_supp_func!(E, nv)
            if val2 == val
                @info "facet confirmed"
                push!(facts_confirmed, fc)
                continue
            end

            w = 1000 * nv + rand(-10:10, length(nv))
            new_vert = elim_vertex!(E, w)

            if !(new_vert in P) # check if new vertex was actually obtained
                @info "new vertex"
                P = convex_hull(P, convex_hull([new_vert]))
                all_confirmed = false
                break
            end
        end
    end

    return P
end

function elim_vertex!(E::ElimData, covec::Vector{Int})
    if E.current_covec != covec
        deform_new_covector!(E, covec)
    end
    return elim_vertex(E)
end

function elim_supp_func!(E::ElimData, covec::Vector{Int})
    if E.current_covec != covec
        deform_new_covector!(E, covec)
    end
    return elim_supp_func(E)
end

function deform_new_covector!(E::ElimData, new_covec::Vector{Int})

    n = size(E.V, 1) - 1
    w_dot_old = (permutedims(vcat(zeros(Int, n), E.current_covec)) * E.A) .- E.shift
    w_dot_new = (permutedims(vcat(zeros(Int, n), new_covec)) * E.A)
    new_shift = min(minimum(w_dot_new) - 1, 0) 
    w_dot_new = w_dot_new .- new_shift

    scl = Int(ceil(maximum(w_dot_new ./ w_dot_old))) + 10
    w_dot_old *= scl
    
    A_start = hcat(vcat(E.A[1:n, :], w_dot_old),
                   vcat(E.A[1:n, :], zeros(Int, 1, size(E.A, 2))))
    A_target = hcat(vcat(E.A[1:n, :], w_dot_new),
                    vcat(E.A[1:n, :], zeros(Int, 1, size(E.A, 2))))

    V_ext = hcat(E.V, E.V)
    wd, p0, p1 = homotopy(A_target, V_ext, A_start, V_ext,
                          E.current_ms, E.current_lift)
    with_logger(NullLogger()) do
        walk_homotopy!(wd, p0, p1)
    end

    sz = size(A_target, 2)
    M_final = MCI(wd.M.V[:, 1:sz], wd.M.A[:, 1:sz])
    E.current_ms = gather_mixed_cells(wd, sz)
    E.current_lift = p1.eps[1:sz]
    E.current_covec = new_covec
    E.shift = new_shift
end

function symbolic_volume(A_mod::Matrix{C},
                         m::MixedCellInds,
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
    
function affine_hull_int(P::Polyhedron{QQFieldElem})
    af = affine_hull(P)
    af_int = AffineHyperplane{QQFieldElem}[]
    for h in af
        g = lcm((denominator).(h.a[1, :])..., denominator(h.b))
        push!(af_int, affine_hyperplane(g * h.a[1, :], g * h.b))
    end
    return af_int
end
              
