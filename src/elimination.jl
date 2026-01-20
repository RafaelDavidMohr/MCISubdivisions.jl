# --- Functions for tropical elimination --- #

function elim_support_func!(E::ElimData,
                            covec::Vector{Int})

    wd = WalkData(E.M_elim, E.current_ms)
    new_lift = get_lifting_vector(E, covec)

    # compute new subdivision
    walk_homotopy!(wd, E.current_lift, new_lift)
    new_ms = gather_mixed_cells(E.M_elim, wd)

    F = prime_field_V(E.M)
    n = size(E.M.V, 1) - 1 # n + 1 input equations
    k = size(E.A, 1) - n - 1
    result = zeros(QQFieldElem, k + 1)
    A_elim = E.A[1:n, :] # last coordinates will be coordinates of eliminant
    for m in new_ms
        proj_mtx = Matrix{QQFieldElem}(undef, k+1, n+k+1)
        for i in 1:k+1
            unit_vec = [j == i ? 1 : 0 for j in 1:k+1]
            lift_unit_vec = get_lifting_vector(E, unit_vec)
            onv_unit_vec = outer_normal_vector(A_elim, m, (QQ).(lift_unit_vec.r))
            proj_mtx[i, :] = vcat(onv_unit_vec, unit_vec)
        end
        _, onv = outer_normal_vector(E.M_elim, m, (QQ).(new_lift.r))
        covec_ext = vcat(onv, covec)
        sp = sortperm(1:size(E.M.A_Fl, 2),
                      rev = true,
                      by = i -> dot(covec_ext, E.M.A_Fl[:, i]))
        V_red = Oscar.echelon_form(matrix(F, E.M.V[:, sp]))
        nz_index = findfirst(!iszero, V_red[end, :])
        result += (vol(m, A_elim) * (proj_mtx * E.A[:, sp[nz_index]]))
    end

    E.current_ms = new_ms
    E.current_lift = new_lift

    return result
end

function get_lifting_vector(E::ElimData, covec::Vector{Int})
    n = size(E.M.V, 1) - 1
    return Lift(vec(permutedims(vcat(zeros(Int, n), covec)) * E.A))
end
