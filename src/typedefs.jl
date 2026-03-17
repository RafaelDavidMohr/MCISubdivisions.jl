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
    r_fl::Float64
    r_P::FqFieldElem
    eps_fl::Float64
    eps_P::FqFieldElem
end

function Base.show(io::IO, a::DualNumber)
    print(io, "$(round(a.r_fl, digits = 3)) + ε ⋅ $(round(a.eps_fl, digits = 3))")
end

function Base.iszero(a::DualNumber)
    return iszero(a.r_P) && iszero(a.eps_P)
end

function Base.:(+)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r_fl + b.r_fl, a.r_P + b.r_P,
                      a.eps_fl + b.eps_fl, a.eps_P + b.eps_P)
end

function Base.:(-)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r_fl - b.r_fl, a.r_P - b.r_P,
                      a.eps_fl - b.eps_fl, a.eps_P - b.eps_P)
end

function Base.:(*)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r_fl*b.r_fl, a.r_P*b.r_P,
                      a.r_fl*b.eps_fl + a.eps_fl*b.r_fl,
                      a.r_P*b.eps_P + a.eps_P*b.r_P)
end

function Base.inv(a::DualNumber)
    iszero(a.r_P) && error("not invertible")
    return DualNumber(a.r_fl^(-1), a.r_P^(-1),
                      -(a.eps_fl*a.r_fl^(-2)), -(a.eps_P*a.r_P^(-2)))
end

function dual_zero(F::FqField)
    return DualNumber(0.0, F(0), 0.0, F(0))
end

function dual_one(F::FqField)
    return DualNumber(1.0, F(1), 0.0, F(0))
end

function lt_dual(a::DualNumber, b::DualNumber)

    if a.r_P != b.r_P
        return a.r_fl < b.r_fl
    elseif a.eps_P != b.eps_P
        return a.eps_fl < b.eps_fl
    else
        return false
    end
end

function prime_field(a::DualNumber)
    return parent(a.r_P)
end

struct DualVector
    R_fl::Vector{Float64}
    R_P::Vector{FqFieldElem}
    E_fl::Vector{Float64}
    E_P::Vector{FqFieldElem}
end

Base.length(v::DualVector) = length(v.R_fl)

function Base.hash(v::DualVector, h::UInt)
    h1 = hash(v.R_P, h)
    return hash(v.E_P, h1)
end

function Base.:(==)(v1::DualVector, v2::DualVector)
    if v1.R_P == v2.R_P
        return v1.E_P == v2.E_P
    end
    return false
end

function prime_field(v::DualVector)
    return parent(first(v.R_P))
end

struct DualMatrix
    R_fl::Matrix{Float64}
    R_P::Matrix{FqFieldElem}
    E_fl::Matrix{Float64}
    E_P::Matrix{FqFieldElem}
end

Base.size(A::DualMatrix, i::Int) = size(A.R_fl, i)

function DualMatrix(F::FqField, A::Matrix{Int}, B::Matrix{Int})
    return DualMatrix((Float64).(A), (F).(A), (Float64).(B), (F).(B))
end

function prime_field(A::DualMatrix)
    return parent(first(A.R_P))
end

# --- Circuit --- #

struct HypVector
    cfs_P::Vector{FqFieldElem}
    cfs_fl::Vector{Float64}
end

Base.hash(v::HypVector, h::UInt) = hash(v.cfs_P, h)
Base.:(==)(v1::HypVector, v2::HypVector) = v1.cfs_P == v2.cfs_P

Base.length(v::HypVector) = length(v.cfs_P)

function prime_field(v::HypVector)
    return parent(first(v.cfs_P))
end

struct Hyperplane{V <: Union{DualVector, HypVector}}
    cfs::V
    nzinds::Vector{Int}

    function Hyperplane(cfs::V) where V <: Union{HypVector, DualVector}
        nzinds = if V <: HypVector
            findall(!iszero, cfs.cfs_P)
        else
            sort(union(findall(!iszero, cfs.R_P), findall(!iszero, cfs.E_P)))
        end
        ni = first(nzinds)
        return new(cfs_P[ni]^(-1) .* cfs_P, cfs_fl[ni]^(-1) .* cfs_fl, nzinds)
   end
end

Base.length(c::Hyperplane) = length(c.cfs_P)

function prime_field(c::Hyperplane)
    return prime_field(c.cfs)
end

# --- MCI --- #

struct SupportMatrix
    A_modP::Matrix{FqFieldElem}
    A_Fl::Matrix{Float64}
end

Base.size(A::SupportMatrix, i::Int) = size(A.A_Fl, i

function prime_field(A::SupportMatrix)
    return parent(first(A.A_modP))
end


struct MCI{M <: Union{DualMatrix, SupportMatrix}}
    V::Matrix{FqFieldElem} # stored over random finite field to speed up computations
    A::M
end

struct RelativeMCI{M <: Union{DualMatrix, SupportMatrix}}
    base::MCI{M}
    A_rel::M
    base_to_rel::Vector{Int}
    rel_to_base::Vector{Int}

    function RelativeMCI(base::MCI{M}, A_rel::M,
                         V_rel::Matrix{FqFieldElem},
                         rel_to_base::Vector{Int}) where {M <: Union{DualMatrix, SupportMatrix}}

        base_to_rel = similar(rel_to_base)
        base_to_rel = zeros(Int, size(base.A, 2))
        for (i, j) in enumerate(rel_to_base)
            base_to_rel[j] = i
        end
        return new(base, A_rel, V_rel, base_to_rel, rel_to_base)
    end
end

function RelativeMCI(M::MCI)
    return RelativeMCI(M, M.A, M.V, indices(M))
end

# TODO: do we need to dispatch over types here?
function MCI(V::Matrix{C}, A::Matrix{Int64}; deform = false) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    F = parent(first(Vp))
    rand_mix = matrix(F, (F).(rand(1:characteristic(F)-1, size(V, 1), size(V, 1))))

    A_modP = reduce_mod_rand_prime(A)
    A_Fl = (Float64).(A)

    supp = if deform
        F = parent(first(A_modP))
        B = rand(-100:100, size(A, 1), size(A, 2))
        DualMatrix(A_Fl, A_modP, (Float64).(B), (F).(B))
    else
        SupportMatrix(A_modP, A_Fl)
    end
    
    return MCI(Matrix(rand_mix * matrix(F, Vp)), supp)
end

function circuits(M::RelativeMCI)
    F = prime_field_V(M)
    matr = matroid_from_matrix_columns(matrix(F, M.V_rel))
    return Oscar.circuits(matr)
end

indices(M::MCI) = collect(1:size(M.A, 2))
indices(M::RelativeMCI) = collect(1:size(M.A_rel, 2))

prime_field_A(M::MCI) = prime_field(M.A)
prime_field_V(M::MCI) = parent(first(M.V))
prime_field_V(M::RelativeMCI) = parent(first(M.V_rel))

function ambient_dim(M::MCI)
    return size(M.A, 1)
end

# --- Lift --- #

struct LiftVector
    r::Vector{Int}
    eps::Vector{Int}
end

function LiftVector(r::Vector{Int})
    l = length(r)
    return LiftVector(r, rand(-10000:10000, l))
end

function test_vector(l::LiftVector; rat=10000)
    return rat*l.r + l.eps
end

# --- WalkData --- # 

struct WalkData{M}
    M::MCI{M}
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

# --- ElimData --- #

mutable struct ElimData
    M::MCI
    M_elim::MCI
    current_ms::Vector{MixedCell}
    current_lift::DualVector
end

function get_elim_start_data(F::Vector{<:MPolyRingElem})
    A, V = get_eci_data(F)
    n = length(F) - 1
    A_elim = A[1:n, :]
    M = MCI(V, A)
    M_elim = MCI(V[1:n, :], A[1:n, :])
    
    init_lift, init_cells = mixed_subdivision(A_elim, M_elim.V)

    return ElimData(M, M_elim, init_cells, init_lift)
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
