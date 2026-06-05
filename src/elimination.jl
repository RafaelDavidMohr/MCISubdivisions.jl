# --- Functions for tropical elimination --- #

function construct_polytope!(E::ElimData)

    n = size(E.V, 1) - 1 # n + 1 input equations
    k = size(E.A, 1) - n - 1
    amb_dim = k + 1 
    @info "computing initial vertices"

    vrt = elim_vertex(E)
    @info "new vertex $(vrt)"
    P = convex_hull([vrt])
    nverts = 1
    while nverts < amb_dim + 1
        @info "dimension $(nverts - 1)"
        af = integer_affine_span(P)
        cfs = rand(-10:10, length(af))
        w = make_smaller(sum(cfs .* af))
        vert = elim_vertex!(E, w)
        fv = first(vertices(P))
        if all(h -> iszero(dot(h, vert - (Int).(fv))), af)
            vert = elim_vertex!(E, -w)
            all(h -> iszero(dot(h, vert - (Int).(fv))), af) && break
        end
        @info "new vertex $(vert)"
        P = convex_hull(P, convex_hull([vert]))
        nverts += 1
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
            # w = make_smaller(10000 * nv + rand(-10:10, length(nv)))
            w = 100 * nv + rand(-3:3, length(nv))
            new_vert = elim_vertex!(E, w)
            if dot(nv, new_vert) == val
                @info "facet confirmed"
                push!(facts_confirmed, fc)
                continue
            end

            if !(new_vert in P) # check if new vertex was actually obtained
                @info "new vertex $(new_vert)"
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
    wd = homotopy(A_target, V_ext, A_start, V_ext,
                  E.current_ms, E.current_lift,
                  Int128.(rand(-LSIZE:LSIZE, size(A_target, 2))))
    
    sz = size(A_target, 2)
    p, ms = try 
        walk_homotopy!(wd)
        wd.p1.eps[1:sz], gather_mixed_cells(wd, sz)
    catch OverflowError
        @info "Overflow, recomputing without deformation"
        pp, mss = mixed_subdivision(A_target, V_ext,
                                    Int128.(rand(-LSIZE:LSIZE, size(A_target, 2))))
        pp.eps, mss
    end

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

function make_smaller(v::Vector{Int})
    mx = maximum((abs).(v))
    if mx > 1000
        nd = ndigits(mx)
        div = 10^(nd - 3)
        return (Int).((round).(v ./ div))
    else
        return v
    end
end

function elim_supp_func!(E::ElimDataDeform, covec::Vector{Int})
    n = size(E.V, 1) - 1
    d = vec(transpose(covec) * E.A[n+1:end, :])
    new_lift = DualVector(d)
    new_ms = deform_subdivision(E.Ap_deform, E.Vp,
                                E.current_ms, E.current_lift,
                                new_lift)

    res = 0
    Ap = E.A[1:n, :]
    Ap_lft = vcat(Ap, permutedims(d))
    for m in new_ms
        w = primitive_normal_vector(Ap_lft, m)
        println(w)
        _, vl = eval_supp_func(Ap_lft, E.V, w, n + 1)
        res += (lifted_volume(m, Ap, d)*vl)
    end

    E.current_lift = new_lift
    E.current_ms = new_ms

    return res
end

function build_projection(E::ElimDataDeform, m::MixedCellInds,
                          covec::Vector{Int},
                          w::Vector{Int})
    n = size(E.V, 1) - 1
    k = size(E.A, 1) - n
    B = rand(-100:100, k, k)
    λ = solve(matrix(QQ, B), QQ.(covec), side = :right)
    Ap = E.A[1:n, :]
    L = permutedims(B) * E.A[n+1:end, :]
    normal_matrix = hcat([normal_space(vcat(Ap, permutedims(L[i, :])), m)[1:n, :]
                          for i in 1:size(L, 1)]...)
    display(normal_matrix)
    μ = solve(matrix(QQ, normal_matrix), QQ.(w), side = :right)
    normal_matrix = normal_matrix * Matrix(diagonal_matrix(μ .// λ))
    return matrix(QQ, normal_matrix) * matrix(QQ, B)^(-1)
end

function elim_vertex!(E::ElimDataDeform, covec::Vector{Int})
    n = size(E.V, 1) - 1
    d = vec(transpose(covec) * E.A[n+1:end, :])
    new_lift = DualVector(d)
    new_ms = deform_subdivision(E.Ap_deform, E.Vp,
                                E.current_ms, E.current_lift,
                                new_lift)

    res = zeros(Int, size(E.A, 1) - n)
    Ap = E.A[1:n, :]
    Ap_lft = vcat(Ap, permutedims(d))
    for m in new_ms
        w = primitive_normal_vector(Ap_lft, m)
        idx, _ = eval_supp_func(Ap_lft, E.V, w, n + 1)
        res += (lifted_volume(m, Ap, d)*E.A[n+1:end, idx])
    end

    E.current_lift = new_lift
    E.current_ms = new_ms

    return res
end

function eval_elim_supp_func(A::Matrix{Int}, V::Matrix{C},
                             covec::Vector{Int},
                             deform::Matrix{Int}) where C

    n = size(V, 1) - 1
    Ap = A[1:n, :]
    Vp = rand(-100:100, n, n + 1) * V
    d = vec(transpose(vcat(zeros(Int, n), covec)) * A)
    Ap_eps = 1000 * Ap + deform
    _, ms = mixed_subdivision(Ap_eps, Vp, d)
    res = 0
    Ap_lft = vcat(Ap, permutedims(d))
    for m in ms
        w = primitive_normal_vector(Ap_lft, m)
        _, vl = eval_supp_func(Ap_lft, V, w, n + 1)
        res += (lifted_volume(m, Ap, d)*vl)
    end
    return res
end
