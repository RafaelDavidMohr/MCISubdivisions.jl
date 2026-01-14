# --- Functions for tropical elimination --- #

function elim_support_func!(E::ElimData,
                            covec::Vector{Int})

    n = size(E.M.V, 1) - 1 # n + 1 input equations
    wd = WalkData(E.M_elim, E.current_ms)
    A_elim = E.A[1:n, :] # last coordinates will be coordinates of eliminant
    new_lift = permutedims(covec) * A_elim 

    # compute new subdivision
    walk_homotopy!(wd, E.current_lift, new_lift)
    new_ms = gather_mixed_cells(wd)

    F = prime_field_V(E.M)
    covec_fl = (Float64).(covec)
    result = 0
    for m in new_ms
        _, onv = outer_normal_vector(E.M_elim, m, (QQ).(new_lift))
        covec_ext = vcat(onv, covec_fl)
        dotprods = permutedims(covec_ext) * E.M.A_fl
        sp = sortperm(dotprods, rev = true)
        V_red = Oscar.echelon_form(matrix(F, E.M.V[:, sp]))
        nz_index = findfirst(!iszero, V_red[end, :])
        result += round(Int, sp[nz_index])*vol(m, A_elim)
    end

    E.current_ms = new_ms
    E.current_lift = new_lift

    return result
end
