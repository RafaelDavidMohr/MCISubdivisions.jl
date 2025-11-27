struct MixedCell
    inds::Vector{Vector{Int}}
end

# Not sure why we have to overload this
function Base.:(==)(m1::MixedCell, m2::MixedCell)
    length(m1.inds) != length(m2.inds) && return false
    for (i, S) in enumerate(m1.inds)
        S != m2.inds[i] && return false
    end
    return true
end

function Base.hash(m::MixedCell, h::UInt)
    return hash(m.inds, h)
end

Base.length(m::MixedCell) = length(m.inds)
function codim(m::MixedCell)
    isempty(m.inds) && return 0
    return sum((s -> length(s) - 1).(m.inds))
end

function Base.show(io::IO, ::MIME"text/plain", m::MixedCell)
    dims = (s -> length(s) - 1).(m.inds)
    print(io, "Mixed cell of dimension $(dims)")
end

struct Circuit
    inds::Vector{Int}
    cfs_modP::Vector{FqFieldElem}
    cfs_fl::Vector{Float64}

    function Circuit(cfs_modP::Vector{FqFieldElem}, cfs_fl::Vector{Float64})
        inds = findall(!iszero, cfs_modP)
        ni = first(inds)
        return new(inds, cfs_modP[ni]^(-1) .* cfs_modP[inds], cfs_fl[ni]^(-1) .* cfs_fl[inds])
    end
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
        @assert size(V, 2) == size(A_modP, 2) "number of coefficients and monomials does not match."
        F = parent(first(V))
        return new(V, A_modP, A_Fl)
    end
end

function MCI(V::Matrix{C}, A::Matrix{Int64}) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    A_modP = reduce_mod_rand_prime(A)
    A_Fl = (Float64).(A)
    
    return MCI(Vp, A_modP, A_Fl)
end

prime_field_A(M::MCI) = parent(first(M.A_modP))
prime_field_V(M::MCI) = parent(first(M.V))

function ambient_dim(M::MCI)
    return size(M.A_modP, 1)
end

struct WalkData
    M::MCI
    walls::Dict{Circuit, Set{Tuple{MixedCell, Int, Bool}}}
end

function WalkData(M::MCI,
                  initial_mixed_cells::Vector{MixedCell})

    walls = Dict{Circuit, Set{Tuple{MixedCell, Int, Bool}}}()
    for m in initial_mixed_cells
        compute_active_walls!(m, M, walls)
    end
    return WalkData(M, walls)
end

mutable struct HomotopyPath
    t_curr::Float64 # current position in path: t_curr * points[1] + (1 - t) * points[2]
    points::Vector{Vector{Float64}}
end

Base.length(path::HomotopyPath) = length(path.points)
Base.getindex(path::HomotopyPath, i::Int) = path.points[i]
is_completed(path::HomotopyPath) = isone(Base.length(path))
