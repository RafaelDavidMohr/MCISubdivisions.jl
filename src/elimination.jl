# --- Functions for tropical elimination --- #

function construct_polytope!(E::ElimData)

    n = size(E.M.V, 1) - 1 # n + 1 input equations
    k = size(E.A, 1) - n - 1
    amb_dim = k + 1 
    @info "computing initial vertices"

    w = rand(-1000:1000, amb_dim)
    P = convex_hull([elim_vertex!(E, w)])
    dm = 0
    while dm < amb_dim 
        @info "dimension $(dm)"
        af = affine_hull(P)
        cfs = rand(-1000:1000, length(af))
        w = sum(cfs .* [intify(h.a[1, :]) for h in af])
        val = sum(cfs .* [intify(h.b) for h in af])
        val2 = elim_supp_func!(E, w)
        w = val == val2 ? -w : w
        vert = elim_vertex!(E, w)
        P = convex_hull(P, convex_hull([vert]))
        dim(P) == dm && break
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

            new_vert = elim_vertex!(E, nv)

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

function elim_vertex!(E::ElimData,
                      covec::Vector{Int})

    covec_d = DualVector(covec)

    covec_mixed_subdivision!(E, covec_d)

    F = prime_field_V(E.M)
    n = size(E.M.V, 1) - 1 # n + 1 input equations
    k = size(E.A, 1) - n - 1
    result = zeros(QQFieldElem, k + 1)
    A_elim = E.A[1:n, :] # last coordinates will be coordinates of eliminant
    for m in E.current_ms
        proj_mtx = Matrix{QQFieldElem}(undef, k+1, n+k+1)
        for i in 1:k+1
            unit_vec = [j == i ? 1 : 0 for j in 1:k+1]
            lift_unit_vec = get_lifting_vector(E, unit_vec)
            onv_unit_vec = outer_normal_vector(A_elim, m, (QQ).(lift_unit_vec))
            proj_mtx[i, :] = vcat(onv_unit_vec, unit_vec)
        end
        _, onv = outer_normal_vector(E.M_elim, m, (QQ).(E.current_lift.r))
        covec_ext = vcat(onv, covec)
        sp = sortperm(1:size(E.M.A_Fl, 2),
                      rev = true,
                      by = i -> dot(covec_ext, E.M.A_Fl[:, i]))
        V_red = Oscar.echelon_form(matrix(F, E.M.V[:, sp]))
        nz_index = findfirst(!iszero, V_red[end, :])
        result += (vol(m, A_elim) * (proj_mtx * E.A[:, sp[nz_index]]))
    end
    return result
end

function elim_supp_func!(E::ElimData,
                         covec::Vector{Int})

    covec_d = DualVector(covec, zeros(Int, length(covec)))

    covec_mixed_subdivision!(E, covec_d)

    F = prime_field_V(E.M)
    n = size(E.M.V, 1) - 1 # n + 1 input equations
    k = size(E.A, 1) - n - 1
    result = QQ(0)
    A_elim = E.A[1:n, :] # last coordinates will be coordinates of eliminant
    for m in E.current_ms
        _, onv = outer_normal_vector(E.M_elim, m, (QQ).(E.current_lift.r))
        covec_ext = (x -> round(Int, x)).(vcat(onv, covec))
        dps = vec(permutedims(covec_ext) * E.A)
        sp = sortperm(dps, rev = true)
        V_red = Oscar.echelon_form(matrix(F, E.M.V[:, sp]), reduced = false)
        nz_index = findfirst(!iszero, V_red[end, :])
        result += (vol(m, A_elim) * dps[sp][nz_index])
    end
    return result
end

function covec_mixed_subdivision!(E::ElimData,
                                  covec::DualVector)

    wd = WalkData(E.M_elim, E.current_ms)
    new_lift = get_lifting_vector(E, covec)
    with_logger(NullLogger()) do
        walk_homotopy!(wd, E.current_lift, new_lift)
    end
    new_ms = gather_mixed_cells(E.M_elim, wd)
    E.current_ms = new_ms
    E.current_lift = new_lift
end

function get_lifting_vector(E::ElimData, covec::Vector{Int})
    n = size(E.M.V, 1) - 1
    return vec(permutedims(vcat(zeros(Int, n), covec)) * E.A)
end

function get_lifting_vector(E::ElimData, covec::DualVector)
    n = size(E.M.V, 1) - 1
    rn = rand(-10000:10000, size(E.A, 2)) 
    return DualVector(get_lifting_vector(E, covec.r),
                      get_lifting_vector(E, covec.eps) + rn)
end

function intify(a::QQFieldElem)
    return Int(denominator(a)*a)
end

function intify(a::Vector{QQFieldElem})
    m = lcm((denominator).(a)...)
    return (Int).(m*a)
end
