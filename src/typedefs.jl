struct MixedCell
    inds::Vector{Vector{Int}}
end

Base.length(m::MixedCell) = length(m.inds)

function Base.show(io::IO, ::MIME"text/plain", m::MixedCell)
    dims = (s -> length(S) - 1).(m.inds)
    print(io, "Mixed cell of dimension $(dims)")
end

struct SparseVec{C}
    cfs::Vector{C}
    inds::Vector{Int}
end

struct Circuit
    inds::Vector{Int}
    cfs_modP::Vector{FqFieldElem}
    cfs_fl::Vector{Float64}
end

function Base.:(==)(c1::Circuit, c2::Circuit)
    return c1.inds == c2.inds && c1.cfs_modP == c2.cfs_modP
end

function Base.hash(c::Circuit, h::UInt)
    return hash(c.inds, hash(c.cfs_modP, h))
end

struct MCI
    V::Matrix{FqFieldElem} # stored over random finite field to speed up computations
    A_modP::Matrix{FqFieldElem}
    A_Fl::Matrix{Float64}

    function MCI(V::Matrix{FqFieldElem}, A_modP::Matrix{FqFieldElem}, A_Fl::Matrix{Float64})
        @assert size(V, 2) == size(A, 2) "number of coefficients and monomials does not match."
        F = base_ring(first(A_modP))
        A_ext_modP = vcat(A_modP, [one(F) for i in 1:1, j in 1:size(A, 2)])
        A_ext_Fl = vcat(A_Fl, ones(Float64, 1, size(A, 2)))
        return new(V, A_ext_modP, A_ext_Fl)
    end
end

function MCI(V::Matrix{QQFieldElem}, A::Matrix{Int64})
    Vp = reduce_mod_rand_prime(V)
    A_modP = reduce_mod_rand_prime(A)
    A_Fl = (Float64).(A)
    
    return MCI(Vp, A_modP, A_Fl)
end

function ambient_dim(M::MCI)
    return size(M.A_modP, 1)
end

struct WalkData
    M::MCI
    walls::Dict{Circuit, Vector{Tuple{Int, MixedCell}}}
end

function WalkData(M::MCI,
                  initial_mixed_cells::Vector{MixedCell})

    walls = Dict{Circuit, Vector{Tuple{Int, MixedCell}}}
    for m in initial_mixed_cells
        compute_active_walls!(m, M, walls)
    end
    return WalkData(M, walls)
end
