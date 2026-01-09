module MCISubdivisions

using Oscar
using LinearAlgebra

include("typedefs.jl")
include("helpers.jl")
include("mixed_subdivs.jl")

export mixed_volume, mixed_subdivision

# --- Main functions --- #

function mixed_volume(F::Vector{<:MPolyRingElem})
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    return mixed_volume(A, V)
end

function mixed_subdivision(F::Vector{<:MPolyRingElem})
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    A, mixed_subdivision(A, V)
end

function mixed_volume(A::Matrix{Int}, V::Matrix{C}) where C
    A_ext, cells = mixed_subdivision(A, V)
    return sum([vol(m, A_ext) for m in cells])
end

function mixed_subdivision(A::Matrix{Int}, V::Matrix{C}) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    F = parent(first(Vp))
    rand_mix = matrix(F, (F).(rand(1:characteristic(F)-1, size(V, 1), size(V, 1))))
    Vp = Matrix(rand_mix * matrix(F, Vp))
    A_ext, wd, p0, p1 = total_degree_homotopy(A, Vp)

    walk_homotopy!(wd, p0, p1)

    A_size = size(A, 2)
    result = Set{MixedCell}()
    for c in keys(wd.walls)
        for (m, act_index, sgn) in wd.walls[c]
            any(S -> any(i -> i > A_size, S), m.inds) && continue
            push!(result, m)
        end
    end

    return A_ext, collect(result)
end

end # module MCISubdivisions
