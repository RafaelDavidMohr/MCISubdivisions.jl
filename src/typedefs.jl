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

# --- Circuit --- #

struct Hyperplane
    cfs_P::Vector{FqFieldElem}
    cfs_fl::Vector{Float64}
    nzinds::Vector{Int}

    function Hyperplane(cfs_P::Vector{FqFieldElem}, cfs_fl::Vector{Float64})
        nzinds = findall(!iszero, cfs_P)
        ni = first(nzinds)
        return new(cfs_P[ni]^(-1) .* cfs_P, cfs_fl[ni]^(-1) .* cfs_fl, nzinds)
   end
end

function Base.length(c::Hyperplane)
    return length(c.cfs_P)
end

function Base.:(==)(c1::Hyperplane, c2::Hyperplane)
    return c1.cfs_P == c2.cfs_P
end

function Base.hash(c::Hyperplane, h::UInt)
    return hash(c.cfs_P, h)
end

function prime_field(c::Hyperplane)
    return parent(first(c.cfs_P))
end

# --- MCI --- #

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

struct RelativeMCI
    base::MCI
    A_rel::Matrix{FqFieldElem}
    V_rel::Matrix{FqFieldElem}
    base_to_rel::Vector{Int}
    rel_to_base::Vector{Int}

    function RelativeMCI(base::MCI, A_rel::Matrix{FqFieldElem},
                         V_rel::Matrix{FqFieldElem},
                         rel_to_base::Vector{Int})

        base_to_rel = similar(rel_to_base)
        base_to_rel = zeros(Int, size(base.A_modP, 2))
        for (i, j) in enumerate(rel_to_base)
            base_to_rel[j] = i
        end
        return new(base, A_rel, V_rel, base_to_rel, rel_to_base)
    end
end

function RelativeMCI(M::MCI)
    return RelativeMCI(M, M.A_modP, M.V, indices(M))
end

function MCI(V::Matrix{C}, A::Matrix{Int64}) where C
    Vp = C <: FqFieldElem ? V : reduce_mod_rand_prime(V)
    A_modP = reduce_mod_rand_prime(A)
    A_Fl = (Float64).(A)
    
    return MCI(Vp, A_modP, A_Fl)
end

function circuits(M::RelativeMCI)
    F = prime_field_V(M)
    matr = matroid_from_matrix_columns(matrix(F, M.V_rel))
    return Oscar.circuits(matr)
end

indices(M::MCI) = collect(1:size(M.A_modP, 2))
indices(M::RelativeMCI) = collect(1:size(M.A_rel, 2))

prime_field_A(M::MCI) = parent(first(M.A_modP))
prime_field_V(M::MCI) = parent(first(M.V))
prime_field_V(M::RelativeMCI) = parent(first(M.V_rel))

function ambient_dim(M::MCI)
    return size(M.A_modP, 1)
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

function test_vector(l::DualVector)
    eps = 10^(-10)
    return l.r + eps .* l.eps
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

# --- ElimData --- #

mutable struct ElimData
    M::MCI
    M_elim::MCI
    A::Matrix{Int}
    current_ms::Vector{MixedCell}
    current_lift::DualVector
end

function get_elim_start_data(F::Vector{<:MPolyRingElem})
    A, V = get_eci_data(F)
    K = base_ring(parent(first(F)))
    rand_mix = matrix(K, (K).(rand(-10000:10000, size(V, 1), size(V, 1))))
    V = Matrix(rand_mix * matrix(K, V))
    n = length(F) - 1
    A_elim = A[1:n, :]
    M = MCI(V, A)
    M_elim = MCI(V[1:n, :], A_elim)
    
    init_lift, init_cells = mixed_subdivision(A_elim, V[1:n, :])

    return ElimData(M, M_elim, A, init_cells, init_lift)
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
