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
        lft = get_lifting_vector(E, covec_d, randomize = false)
        onv = outer_normal_vector(A_elim, m, (QQ).(test_vector(lft)))
        covec_ext = vcat(onv, test_vector(covec_d))
        sp = sortperm(1:size(E.A, 2),
                      rev = true,
                      by = i -> dot(covec_ext, E.A[:, i]))
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
    for m in E.current_ms
        Al = vcat(E.A[1:n, :], permutedims(E.current_lift.r))
        println("cell indices $(m.inds)")
        ns = normal_space(Al, m.inds)
        println("normal space: ")
        display(ns)
        if size(ns, 2) > 1
            println("---")
            continue
        end
        troproot = (Int).(lcm((denominator).(ns[:, 1]))*ns[:, 1])
        if !iszero(last(troproot))
            ns *= last(troproot)^(-1)
        end
        dps = vec(permutedims(troproot) * Al)
        println("dot products: $(dps)")
        _, i = findmax(dps)
        if !(i in first(m.inds))
            if !iszero(last(troproot))
                println("---")
                continue
            end
            troproot = -troproot
            dps = -dps
            _, i = findmax(dps)
            if !(i in first(m.inds))
                println("---")
                continue
            end
        end
        sp = sortperm(dps, rev = true)
        println("sorted indss: $(sp)")
        V_red = Oscar.echelon_form(matrix(F, E.M.V[:, sp]), reduced = false)
        println("rank $(rank(matrix(F, E.M.V[:, first(m.inds)])))")
        nz_index = findfirst(!iszero, V_red[end, :])
        result += (vol(m.inds, E.A) * dps[sp][nz_index])
        println("---")
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

function get_lifting_vector(E::ElimData, covec::DualVector; randomize=true)
    n = size(E.M.V, 1) - 1
    rn = rand(-10000:10000, size(E.A, 2)) 
    epsp = get_lifting_vector(E, covec.eps)
    return DualVector(get_lifting_vector(E, covec.r),
                      randomize ? epsp + rn : epsp)
end

function intify(a::QQFieldElem)
    return Int(denominator(a)*a)
end

function intify(a::Vector{QQFieldElem})
    m = lcm((denominator).(a)...)
    return (Int).(m*a)
end
