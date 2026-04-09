module MCISubdivisions

using Oscar
using LinearAlgebra
using Logging

include("typedefs.jl")
include("helpers.jl")
include("mixed_subdivs.jl")
include("elimination.jl")

export mixed_volume, mixed_subdivision, get_eliminant_polytope, get_A_disc_equations

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
    return A, mixed_subdivision(A, V)[2]
end

function mixed_volume(A::Matrix{Int}, V::Matrix{C}) where C
    _, cells = mixed_subdivision(A, V)
    return sum([vol(m, A) for m in cells])
end

function mixed_subdivision(A::Matrix{Int}, V::Matrix{C}) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    wd, p0, p1 = starting_system(A, Vp)

    walk_homotopy!(wd, p0, p1)

    A_size = size(A, 2)
    M_final = MCI(wd.M.V[:, 1:A_size], wd.M.A[:, 1:A_size])
    return DualVector(p1.r[1:A_size], p1.eps[1:A_size]), gather_mixed_cells(M_final, wd, A_size)
end

function get_eliminant_polytope(F::Vector{<:MPolyRingElem})
    A, V = get_eci_data(F)
    Vp = eltype(V) <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    @info "setting up initial data"
    E = ElimData(A, Vp)
    @info "done"
    return construct_polytope!(E)
end

end # module MCISubdivisions
