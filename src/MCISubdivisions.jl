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

function mixed_volume(F::Vector{<:MPolyRingElem}; lift_type = Int, arithmetic = :ff)
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    if arithmetic == :ff
        Vp = typeof(first(F)) <: FqMPolyRingElem ? V : reduce_mod_rand_prime(V)
        return mixed_volume(A, Vp, lift_type = lift_type)
    elseif arithmetic == :float
        return mixed_volume(A, V, lift_type = lift_type)
    else
        error("unrecognized keyword value $(arithmetic)")
    end
end

function mixed_subdivision(F::Vector{<:MPolyRingElem};
                           lift_type = Int, arithmetic = :ff)
    
    R = parent(first(F))
    @assert ngens(R) == length(F) "Input system not square"
    A, V = get_eci_data(F)
    d = lift_type.(rand(-LSIZE:LSIZE, size(A, 2)))
    if arithmetic == :ff
        Vp = typeof(first(F)) <: FqMPolyRingElem ? V : reduce_mod_rand_prime(V)
        return A, mixed_subdivision(A, Vp, d)[2]
    elseif arithmetic == :float
        return A, mixed_subdivision(A, V, d)[2]
    else
        error("unrecognized keyword value $(arithmetic)")
    end
end

function mixed_volume(A::Matrix{Int}, V::Matrix{C}; lift_type = Int) where C

    d = lift_type.(rand(-LSIZE:LSIZE, size(A, 2)))
    _, cells = mixed_subdivision(A, V, d)
    return sum([vol(m, A) for m in cells])
end

function mixed_subdivision(A::Matrix{Int}, V::Matrix{C},
                           d::Vector{T}) where {C, T}
    
    @assert C <: Int || C <: FqFieldElem "Only integer of finite field coefficients supported"

    wd = starting_system(A, V, d)

    walk_homotopy!(wd)

    p1 = wd.p1
    A_size = size(A, 2)
    return p1[1:A_size], gather_mixed_cells(wd, A_size)
end

function real_root_count(A::Matrix{Int}, V::Matrix{QQFieldElem}, ms::Vector{MixedCellInds})
    return sum([real_root_count(m, A, V) for m in ms])
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
