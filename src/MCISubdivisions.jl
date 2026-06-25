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

export mixed_volume, mixed_subdivision, real_root_count, get_eliminant_polytope, get_A_disc_equations

# --- Main functions --- #

"""
    mixed_volume(F::Vector{<:MPolyRingElem})

Compute the mixed volume of the engineered complete intersection
defined by (V, A) where A is the support and V the coefficient matrix
of `F`.
"""
function mixed_volume(F::Vector{<:MPolyRingElem})
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    Vp = typeof(first(F)) <: FqMPolyRingElem ? V : reduce_mod_rand_prime(V)
    return mixed_volume(A, Vp)
end

"""
    mixed_subdivision(F::Vector{<:MPolyRingElem})

Compute a mixed subdivision at a randomly generated height vector of
the engineered complete intersection defined by (V, A) where A is the
support and V the coefficient matrix of `F`.
"""
function mixed_subdivision(F::Vector{<:MPolyRingElem})
    
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    d = rand(-LSIZE:LSIZE, size(A, 2))
    Vp = reduce_mod_rand_prime(V)
    return A, mixed_subdivision(A, Vp, d)[2]
end

"""
    mixed_volume(A::Matrix{Int}, V::Matrix{C}) where C

Compute the mixed volume of the engineered complete intersection
defined by `V` and `A`.
"""
function mixed_volume(A::Matrix{Int}, V::Matrix{C}) where C

    d = rand(-LSIZE:LSIZE, size(A, 2))
    _, cells = mixed_subdivision(A, V, d)
    return sum([vol(m, A) for m in cells])
end

"""
    mixed_subdivision(A::Matrix{Int}, V::Matrix{C},
                      d::Vector{Int} = rand(-LSIZE:LSIZE, size(A, 2)),
                      d_start::Vector{Int} = rand(-LSIZE:LSIZE, size(A, 2))) where C

Compute the mixed subdivision of the engineered complete intersection
defined by `V` and `A` at the height vector `d` starting from the
regular subdivision of `A` defined by `d_start`. This may fail if
`d_start` is not sufficiently generic.
"""
function mixed_subdivision(A::Matrix{Int}, V::Matrix{C},
                           d::Vector{Int} = rand(-LSIZE:LSIZE, size(A, 2)),
                           d_start::Vector{Int} = rand(-LSIZE:LSIZE, size(A, 2))) where C
    
    @assert C <: Int || C <: FqFieldElem "Only integer or finite field coefficients supported"

    wd = starting_system(A, V, d, d_start)

    walk_homotopy!(wd)

    p1 = wd.p1
    A_size = size(A, 2)
    return p1[1:A_size], gather_mixed_cells(wd, A_size)
end

"""
    real_root_count(A::Matrix{Int}, V::Matrix{QQFieldElem},
                    ms::Vector{MixedCellInds})

Compute the real root count of `A`, `V` and the mixed subdivision given by `ms`.
This is defined to be the sum of the real root counts of the polynomial systems
corresponding to the leading parts of `A` and `V` indicated by `ms`.
"""
function real_root_count(A::Matrix{Int}, V::Matrix{QQFieldElem},
                         ms::Vector{MixedCellInds})

    return sum([real_root_count(m, A, V) for m in ms])
end

"""
    get_eliminant_polytope(F::Vector{<:MPolyRingElem})

Compute the Newton polytope of the eliminant of the engineered
complete intersection (A, V) where `A` is the support and `V` the
coefficient matrix of `F`. The last k variables are eliminated where
`F` consists of k + 1 equations.
"""
function get_eliminant_polytope(F::Vector{<:MPolyRingElem})
    A, V = get_eci_data(F)
    @info "setting up initial data"
    E = ElimData(A, V)
    @info "done"
    return construct_polytope!(E)
end

end # module MCISubdivisions
