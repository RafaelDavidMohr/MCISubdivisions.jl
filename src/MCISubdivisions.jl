module MCISubdivisions

using Oscar
using LinearAlgebra
using Logging
using SparseArrays
using DataStructures

include("typedefs.jl")
include("helpers.jl")
include("mixed_subdivs.jl")
include("elimination.jl")

export mixed_volume, mixed_subdivision, get_eliminant_polytope, get_A_disc_equations

# --- Main functions --- #

function mixed_volume(F::Vector{<:MPolyRingElem}; lift_type = Int)
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    return mixed_volume(A, V, lift_type = lift_type)
end

function mixed_subdivision(F::Vector{<:MPolyRingElem}; lift_type = Int)
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    d = lift_type.(rand(-LSIZE:LSIZE, size(A, 2)))
    return A, mixed_subdivision(A, V, d)[2]
end

function mixed_volume(A::Matrix{Int}, V::Matrix{C}; lift_type = Int) where C

    d = lift_type.(rand(-LSIZE:LSIZE, size(A, 2)))
    _, cells = mixed_subdivision(A, V, d)
    return sum([vol(m, A) for m in cells])
end

function mixed_subdivision(A::Matrix{Int}, V::Matrix{C},
                           d::Vector{T}) where {C, T}
    
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    wd = starting_system(A, Vp, d)

    walk_homotopy!(wd)

    p0, p1 = wd.p0, wd.p1
    A_size = size(A, 2)
    return DualVector(p1.r[1:A_size], p1.eps[1:A_size]), gather_mixed_cells(wd, A_size)
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
