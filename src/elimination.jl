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
        dps, troproot = get_tropical_root(E.A[1:n, :], E.current_lift.r, m.inds)
        isnothing(troproot) && continue
        proj_mtx = Matrix{QQFieldElem}(undef, k+1, n+k+1)
        for i in 1:k+1
            unit_vec = [j == i ? 1 : 0 for j in 1:k+1]
            lift_unit_vec = get_lifting_vector(E, unit_vec)
            Al = vcat(E.A[1:n, :], permutedims(lift_unit_vec))
            ns = normal_space(Al, m.inds)
            size(ns, 2) > 1 && error("don't know what to do in this case")
            onv_unit_vec = ns[:, 1]
            if iszero(last(onv_unit_vec))
                proj_mtx[i, :] = vcat(onv_unit_vec[1:end-1], zeros(Int, k+1))
            else
                onv_unit_vec *= last(onv_unit_vec)^(-1)
                proj_mtx[i, :] = vcat(onv_unit_vec[1:end-1], unit_vec)
            end
        end
        troproot_ext = if iszero(last(troproot))
            vcat(troproot[1:end-1], zeros(Int, length(covec)))
        else
            vcat(troproot[1:end-1], covec)
        end
        sol = solve(matrix(QQ, transpose(proj_mtx)),
                    (QQ).(troproot_ext), side = :right)
        for i in 1:length(sol)
            if sol[i] != covec[i]
                proj_mtx[i, :] = -proj_mtx[i, :]
            end
        end

        sp = sortperm(dps, rev = true)
        V_red = Oscar.echelon_form(matrix(F, E.M.V[:, sp]), reduced = false)
        nz_index = findfirst(!iszero, V_red[end, :])
        result += (vol(m.inds, E.A) * (proj_mtx * E.A[:, sp[nz_index]]))
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
        dps, troproot = get_tropical_root(E.A[1:n, :], E.current_lift.r, m.inds)
        if isnothing(troproot)
            println("---")
            continue
        end
        println("is non-trivial canc")
        sp = sortperm(dps, rev = true)
        println("sorted :   $(sp)")
        println("dotps  :   $(dps)")
        V_red = Oscar.echelon_form(matrix(F, E.M.V[:, sp]), reduced = false)
        nz_index = findfirst(!iszero, V_red[end, :])
        println("volume :   $(vol(m.inds, E.A))")
        println("las ind:   $(sp[nz_index])")
        result += (vol(m.inds, E.A) * dps[sp][nz_index])
        println("---")
    end
    return result
end

function get_tropical_root(A::Matrix{Int64}, lift::Vector{Int},
                           m_inds::Vector{Vector{Int}})

    println("indices:   $(m_inds)")
    Al = vcat(A, permutedims(lift))
    ns = normal_space(Al, m_inds)
    if size(ns, 2) > 1
        println("dim too large")
        return nothing, nothing
    end

    troproot = lcm((denominator).(ns[:, 1]))*ns[:, 1]
    println("trop root: $(troproot)")
    # if !iszero(last(troproot))
    #     troproot *= last(troproot)^(-1)
    # end

    dps = vec(permutedims(troproot) * Al)
    _, i = findmax(dps)
    if !(i in first(m_inds))
        if !iszero(last(troproot))
            println("case 1")
            println("dotps  :   $(dps)")
            return nothing, nothing
        end
        troproot = -troproot
        dps = -dps
        _, i = findmax(dps)
        if !(i in first(m_inds))
            println("case 2")
            println("dotps  :   $(dps)")
            return nothing, nothing
        end
    end

    return dps, (Int).(troproot)
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
