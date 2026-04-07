# --- Mixed Cell --- #

struct MixedCell
    inds::Vector{Vector{Int}}
    loc_inds::Vector{Vector{Int}}
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

# --- Dual Number --- #

struct DualNumber
    r::Float64
    eps::Float64
end

function Base.show(io::IO, a::DualNumber)
    print(io, "$(round(a.r, digits = 2)) + ε ⋅ $(round(a.eps, digits = 2))")
end

function Base.iszero(a::DualNumber)
    return iszero(a.r) && iszero(a.eps)
end

function Base.zero(::Type{DualNumber})
    return DualNumber(0.0, 0.0)
end

function Base.one(::Type{DualNumber})
    return DualNumber(1.0, 0.0)
end

function Base.:(+)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r + b.r, a.eps + b.eps)
end

function Base.:(-)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r - b.r, a.eps - b.eps)
end

function Base.:(*)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r*b.r, a.r_P*b.r_P, a.r*b.eps + a.eps*b.r)
end

function Base.inv(a::DualNumber)
    iszero(a.r) && error("not invertible")
    return DualNumber(a.r^(-1), -(a.eps*a.r^(-2)))
end

function lt_dual(a::DualNumber, b::DualNumber)

    if a.r != b.r
        return a.r < b.r
    elseif a.eps != b.eps
        return a.eps < b.eps
    else
        return false
    end
end

# --- Circuit --- #

struct Hyperplane
    cfs::Vector{Int}
    nzinds::Vector{Int}

    function Hyperplane(cfs::Vector{Int})
        nzinds = findall(!iszero, cfs)
        ni = first(nzinds)
        sb = signbit(cfs[ni])
        sb ? return new(-cfs, nzinds) : return(cfs, nzinds)
   end
end

function Base.length(c::Hyperplane)
    return length(c.cfs)
end

function Base.:(==)(c1::Hyperplane, c2::Hyperplane)
    return c1.cfs == c2.cfs
end

function Base.hash(c::Hyperplane, h::UInt)
    return hash(c.cfs, h)
end

# --- MCI --- #

struct MCI
    V::Matrix{FqFieldElem} # stored over random finite field to speed up computations
    A::Matrix{Int}

    function MCI(V::Matrix{FqFieldElem}, A::Matrix{Int})
        @assert size(V, 2) == size(A_modP, 2) "number of coefficients and monomials does not match."
        F = parent(first(V))
        R = echelon_form(matrix(F, V))
        return new(Matrix(R), A)
    end
end

struct RelativeMCI
    base::MCI
    V_rel::Matrix{FqFieldElem}
    base_to_rel::Vector{Int}
    rel_to_base::Vector{Int}

    function RelativeMCI(base::MCI, V_rel::Matrix{FqFieldElem},
                         rel_to_base::Vector{Int})

        base_to_rel = similar(rel_to_base)
        base_to_rel = zeros(Int, size(base.A, 2))
        for (i, j) in enumerate(rel_to_base)
            base_to_rel[j] = i
        end
        return new(base, V_rel, base_to_rel, rel_to_base)
    end
end

function RelativeMCI(M::MCI)
    return RelativeMCI(M, M.V, indices(M))
end

function MCI(V::Matrix{C}, A::Matrix{Int64}) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    F = parent(first(Vp))
    rand_mix = matrix(F, (F).(rand(1:characteristic(F)-1, size(V, 1), size(V, 1))))
    
    return MCI(Matrix(rand_mix * matrix(F, Vp)), A)
end

function circuits(M::RelativeMCI)
    F = prime_field_V(M)
    matr = matroid_from_matrix_columns(matrix(F, M.V_rel))
    return Oscar.circuits(matr)
end

indices(M::MCI) = collect(1:size(M.A, 2))
indices(M::RelativeMCI) = collect(1:size(M.V_rel, 2))

prime_field_A(M::MCI) = parent(first(M.A_modP))
prime_field_V(M::MCI) = parent(first(M.V))
prime_field_V(M::RelativeMCI) = parent(first(M.V_rel))

function ambient_dim(M::MCI)
    return size(M.A, 1)
end

# --- Lift --- #

struct DualVector
    r::Vector{Int}
    eps::Vector{Int}
end

function DualVector(r::Vector{Int})
    l = length(r)
    return DualVector(r, rand(-10000:10000, l))
end

function test_vector(l::DualVector; rat=10000)
    return rat*l.r + l.eps
end

# --- WalkData --- # 

struct WalkData
    M::MCI
    walls::Dict{Hyperplane, Set{Tuple{MixedCell, Int, Bool}}}
    nocross::Set{Hyperplane}
end

function WalkData(M::MCI,
                  initial_mixed_cells::Vector{MixedCell})

    walls = Dict{Hyperplane, Set{Tuple{MixedCell, Int, Bool}}}()
    for m in initial_mixed_cells
        compute_active_walls!(m, M, walls)
    end
    return WalkData(M, walls, Set{Hyperplane}())
end

# --- RRCounter (for convenience) --- #

struct RRCounter
    A::Matrix{Int}
    V::Matrix{QQFieldElem}
    cnt::Dict{MixedCell, Int}
    target_rr_count::Int

    function RRCounter(A::Matrix{Int}, V::Matrix{QQFieldElem}, target_rr_count::Int)
        return new(A, V, Dict{MixedCell, Int}(), target_rr_count)
    end
end

function empty_rr_counter()
    return RRCounter(Matrix{Int}(undef, 0, 0), Matrix{QQFieldElem}(undef, 0, 0), -1)
end

function rr_count(rr_counter::RRCounter)
    return sum([rr_counter.cnt[m] for m in keys(rr_counter.cnt)])
end

function add_mixed_cell!(rr_counter::RRCounter, m::MixedCell)
    rr_counter.cnt[m] = real_root_count(m, rr_counter.A, rr_counter.V)
end
