# --- Functions for tropical elimination --- #

function construct_polytope!(E::Union{ElimData, ElimDataDeform})

    n = size(E.V, 1) - 1 # n + 1 input equations
    k = size(E.A, 1) - n - 1
    amb_dim = k + 1 
    @info "computing initial vertices"

    vrt = elim_vertex(E)
    @info "new vertex $vrt"
    P = convex_hull([vrt])
    nverts = 1
    while nverts < amb_dim + 1
        @info "$(nverts) vertices computed"
        ns = normal_space(P)
        w = int_approximate(ns * rand(-10:10, size(ns, 2))) + rand(-10:10, size(ns, 1))
        vert = elim_vertex!(E, w)
        if vert in P
            vert = elim_vertex!(E, -w)
            vert in P && break
        end
        @info "new vertex $vert"
        P = convex_hull(P, convex_hull([vert]))
        nverts += 1
    end
    @info "done"

    facts_confirmed = Set{AffineHalfspace{QQFieldElem}}()
    vs = Set([Int.(v) for v in vertices(P)])

    all_confirmed = false
    while !all_confirmed
        @info "recomputing facets"
        facts = facets(P)
        @info "done, $(length(facts)) facets to check"
        all_confirmed = true

        for fc in facts
            fc in facts_confirmed && continue

            nv = ZZ.(fc.a[1,:])
            val = fc.b
            w = shrink(nv)
            @info "trying to compute new vertex"
            new_vert = elim_vertex!(E, w)
            if dot(nv, new_vert) == val
                @info "facet confirmed"
                push!(facts_confirmed, fc)
                continue
            else
                all_confirmed = false
            end

            if !(new_vert in vs) # check if new vertex was actually obtained
                @info "new vertex $(new_vert)"
                push!(vs, new_vert)
            end
        end

        P = convex_hull(collect(vs), non_redundant = true)
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

function new_covector!(E::ElimData, new_covec::Vector{Int})
    n = size(E.V, 1) - 1
    w_dot_new = (permutedims(vcat(zeros(Int, n), new_covec)) * E.A)
    new_shift = min(minimum(w_dot_new) - 1, 0) 
    w_dot_new = w_dot_new .- new_shift
    
    A_target = hcat(vcat(E.A[1:n, :], w_dot_new),
                    vcat(E.A[1:n, :], zeros(Int, 1, size(E.A, 2))))

    V_ext = hcat(E.V, E.V)
    
    pp, ms = mixed_subdivision(A_target, V_ext,
                              rand(-LSIZE:LSIZE, size(A_target, 2)))
    p = [pi.eps for pi in pp]

    E.current_ms = ms
    E.current_lift = p
    E.current_covec = new_covec
    E.shift = new_shift
end

# this currently does not work sometimes
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
    
    sz = size(A_target, 2)
    p, ms = try
        wd = homotopy(A_target, V_ext, A_start, V_ext,
                      E.current_ms, E.current_lift,
                      rand(-LSIZE:LSIZE, size(A_target, 2)))
        walk_homotopy!(wd)
        [pi.eps for pi in wd.p1[1:sz]], gather_mixed_cells(wd, sz)
    catch RoundingError
        @info "Rounding error, recomputing without deformation"
        pp, mss = mixed_subdivision(A_target, V_ext,
                                    rand(-LSIZE:LSIZE, size(A_target, 2)))
        [pi.eps for pi in pp], mss
    end

    any(m -> !certify_mixed_cell(A_target, V_ext, m, p), ms) && error("Error in mixed subdivision computation, possibly due to rounding. Try to relaunch the computation.")

    E.current_ms = ms
    E.current_lift = p
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

function integer_affine_span(P::Polyhedron{QQFieldElem})
    vs = vertices(P)
    if isone(length(vs))
        n = Oscar.ambient_dim(P)
        return [[j == i ? 1 : 0 for j in 1:n] for i in 1:n]
    end
    mat = (Int).(vcat([transpose(vs[1] - v) for v in vs[2:end]]...))
    k = kernel(matrix(ZZ, mat), side = :right)
    return filter(!isempty, [(Int).(k[:, i]) for i in 1:size(k, 2)])
end

function normal_space(P::Polyhedron{QQFieldElem})
    vs = vertices(P)
    if isone(length(vs))
        n = Oscar.ambient_dim(P)
        return Float64.(id_matrix(n))
    end
    mat = (Int).(vcat([transpose(vs[1] - v) for v in vs[2:end]]...))
    return nullspace(mat)
end

function compute_roughly_orthogonal_vector(P::Polyhedron{QQFieldElem})
    ns = normal_space(P)
    return Int.(round.(ns * rand(-1000:1000, size(ns, 2))))
end

function elim_supp_func!(E::ElimDataDeform, covec::Vector{Int}; need_deform = true)
    n = size(E.V, 1) - 1
    d = vec(transpose(covec) * E.A[n+1:end, :])

    if need_deform
        new_lift = DualVector(d)
        new_ms = deform_subdivision(E.Ap_deform, E.Vp,
                                    E.current_ms, E.current_lift,
                                    new_lift)
        E.current_lift = new_lift
        E.current_ms = new_ms
    end

    res = 0
    Ap = E.A[1:n, :]
    Ap_lft = vcat(Ap, permutedims(d))
    for m in E.current_ms
        # vl = volume_in_affine_span(m, E.A)
        vl = lifted_volume(m, Ap, d)
        # if iszero(vl)
        #     @warn "warning: zero volume for covec $(covec)"
        #     continue
        # end
        w = primitive_normal_vector(Ap_lft, m)
        _, evl = eval_supp_func(Ap_lft, E.V, w, n + 1)
        res += (vl*evl)
    end

    return res
end

function elim_vertex!(E::ElimDataDeform, covec::Vector{Int})

    n = size(E.V, 1) - 1
    k = size(E.A, 1) - n

    elim_supp_func!(E, covec)

    mat = vcat([permutedims(1000 * covec) for _ in 1:7]...)
    mat += rand(-10:10, k, k)
    rhs = [elim_supp_func!(E, mat[i, :], need_deform=false) for i in 1:size(mat, 1)]

    return Int.(round.(mat \ rhs))
end

function elim_vertex(E::ElimDataDeform)
    n = size(E.V, 1) - 1
    k = size(E.A, 1) - n
    w = rand(-10:10, k)
    return elim_vertex!(E, w)
end

function int_approximate(v::Vector{Float64}, precision = 2)
    return Int.(round.(10^precision * v))
end

function shrink(v::Vector{ZZRingElem})
    vt = BigInt.(v)
    nd = ndigits(maximum(abs.(vt)))
    if nd > 4
        vnew = vt ./ 10^(nd - 4)
        return Int.(round.(vnew))
    else
        return Int.(vt)
    end
end
