# --- Functions for tropical elimination --- #

function elim_support_func(M::MCI,
                           M_elim::MCI,
                           A_int::Matrix{Int},
                           covec::Vector{Int},
                           previous_ms::Vector{MixedCell},
                           previous_lift::Vector{Int})

    n = size(M.V, 1) - 1 # n + 1 input equations
    wd = WalkData(M_elim, previous_ms)
    A_elim = A_int[:, 1:n] # last coordinates will be coordinates of eliminant
    new_lift = permutedims(covec) * A_elim 

    # compute new subdivision
    walk_homotopy!(wd, previous_lift, new_lift)
    new_ms = gather_mixed_cells(wd)

    F = prime_field_V(M)
    covec_fl = (Float64).(covec)
    result = 0
    for m in new_ms
        _, onv = outer_normal_vector(M_elim, m, (QQ).(new_lift))
        covec_ext = vcat(onv, covec_fl)
        dotprods = permutedims(covec_ext) * M.A_fl
        sp = sortperm(dotprods, rev = true)
        V_red = Oscar.echelon_form(matrix(F, M.V[:, sp]))
        nz_index = findfirst(!iszero, V_red[end, :])
        result += round(Int, sp[nz_index])*vol(m, A_elim)
    end

    return result
end
