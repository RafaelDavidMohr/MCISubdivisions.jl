# --- MinExtractor --- #

struct MinExtractor{T}
    heap::BinaryMinHeap{T}
    seen::Set{T}
    
    MinExtractor{T}() where T = new(BinaryMinHeap{T}(), Set{T}())
end

Base.isempty(me::MinExtractor) = isempty(me.seen)

function Base.push!(me::MinExtractor{T}, x::T) where T
    x in me.seen && return
    push!(me.seen, x)
    push!(me.heap, x)
    return
end

function extract_min!(me::MinExtractor{T}) where T
    x = pop!(me.heap)
    delete!(me.seen, x)
    return x
end

# --- Mixed Cell --- #

const MixedCellInds = Vector{Vector{Int}}

struct MixedCell
    inds::MixedCellInds
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

const PRIME = 2^31 - 1

struct FPNum
    x::Int
end

Base.zero(::Type{FPNum}) = FPNum(0)
Base.zero(::FPNum)= zero(FPNum)
Base.iszero(a::FPNum) = iszero(a.x)
Base.one(::Type{FPNum}) = FPNum(1)
Base.one(::FPNum) = one(FPNum)

Base.:(+)(a::FPNum, b::FPNum) = FPNum(mod(a.x + b.x, PRIME))
Base.:(-)(a::FPNum, b::FPNum) = FPNum(mod(a.x + (PRIME - b.x), PRIME))
Base.:(-)(a::FPNum) = FPNum(PRIME - a.x)
Base.:(*)(a::FPNum, b::FPNum) = FPNum(mod(a.x * b.x, PRIME))
Base.:(*)(a::Integer, b::FPNum) = FPNum(mod(a * b.x, PRIME))
Base.:(*)(a::FPNum, b::Integer) = b * a
Base.inv(a::FPNum) = FPNum(invmod(a.x, PRIME))

function Base.:(*)(a::FPNum, v::SparseVector{FPNum, Int})
    res = similar(v)
    for i in findall(!iszero, v)
        res[i] = a * v[i]
    end
    return res
end

const INV = Union{AbstractFloat, Rational, FPNum}

struct DualNumber{T}
    r::T
    eps::T
end

function Base.show(io::IO, a::DualNumber{<:Real})
    print(io, "$(round(a.r, digits = 2)) + ε ⋅ $(round(a.eps, digits = 2))")
end

function Base.iszero(a::DualNumber)
    return iszero(a.r) && iszero(a.eps)
end

function Base.zero(::Type{DualNumber{T}}) where T
    return DualNumber{T}(zero(T), zero(T))
end
Base.zero(::DualNumber{T}) where T = zero(DualNumber{T})

function Base.one(::Type{DualNumber{T}}) where T
    return DualNumber{T}(one(T), zero(T))
end
Base.one(::DualNumber{T}) where T = one(DualNumber{T})

function Base.:(+)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r + b.r, a.eps + b.eps)
end

function Base.:(-)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r - b.r, a.eps - b.eps)
end
function Base.:(-)(a::DualNumber)
    return DualNumber(-a.r, -a.eps)
end

function Base.:(*)(a::DualNumber, b::DualNumber)
    return DualNumber(a.r*b.r, a.r*b.eps + a.eps*b.r)
end
Base.:(*)(a::Union{Real, FPNum}, b::DualNumber) = DualNumber(a*b.r, a*b.eps)
Base.:(*)(b::DualNumber, a::Union{Real, FPNum}) = a * b

isinvertible(a::DualNumber{<:INV}) = !iszero(a.r)
function Base.inv(a::DualNumber{<:INV})
    iszero(a.r) && error("$(a.r), $(a.eps) not invertible")
    arinv = Base.inv(a.r)
    return DualNumber(arinv, -(a.eps*arinv*arinv))
end
function Base.:(/)(a::DualNumber{<:INV}, b::DualNumber{<:INV})
    return a * Base.inv(b)
end

# needed for other linear algebra operations
Base.transpose(a::DualNumber) = a

# needed to compute dot products of vectors
LinearAlgebra.dot(a::DualNumber, b::DualNumber) = a * b
LinearAlgebra.dot(a::Union{Real, FPNum}, b::DualNumber) = a * b
LinearAlgebra.dot(a::DualNumber, b::Union{Real, FPNum}) = b * a


function Base.isless(a::DualNumber{<:Real}, b::DualNumber{<:Real})

    if a.r != b.r
        return a.r < b.r
    elseif a.eps != b.eps
        return a.eps < b.eps
    else
        return false
    end
end

struct OrderDual
    x::DualNumber{Float64}
    xP::DualNumber{FPNum}
end

Base.:(==)(a::OrderDual, b::OrderDual) = a.xP == b.xP

Base.zero(::Type{OrderDual}) = OrderDual(zero(DualNumber{Float64}), zero(DualNumber{FPNum}))
Base.one(::Type{OrderDual}) = OrderDual(one(DualNumber{Float64}), one(DualNumber{FPNum}))

function Base.isless(a::OrderDual, b::OrderDual)
    a.xP == b.xP && return false
    if a.xP.r != b.xP.r
        return a.x.r < b.x.r
    elseif a.xP.eps != b.xP.eps
        return a.x.eps < b.x.eps
    else
        return false
    end
end

# --- Lift --- #

const LSIZE = 50000

const DualVector{T} = Vector{DualNumber{T}} 

function DualVector(r::Vector{<:Integer})
    return [DualNumber(ri, rand(-LSIZE:LSIZE)) for ri in r]
end

# --- Hyperplane --- #

struct Hyperplane
    cfs::SparseVector{Float64, Int}
    cfs_modP::SparseVector{FPNum, Int}
    dot0::DualNumber{Float64}
    dot1::DualNumber{Float64}
    dot0_modP::DualNumber{FPNum}
    dot1_modP::DualNumber{FPNum}
    act_index::Int
    sgn::Bool
    exchange_index::Int
end

function Base.:(==)(c1::Hyperplane, c2::Hyperplane)
    c1.cfs_modP == c2.cfs_modP
end

is_exchange(c::Hyperplane) = !iszero(c.exchange_index)


# --- MCI --- #

struct MCI{C}
    V::Matrix{C}
    A::Matrix{Int}
    rank_cache::Dict{Vector{Int}, Int}

    function MCI{C}(V::Matrix{C}, A::Matrix{Int}) where C
        @assert size(V, 2) == size(A, 2) "number of coefficients and monomials does not match."
        return new(V, A, Dict{Vector{Int}, Int}())
    end
end

MCI(V::Matrix{C}, A::Matrix{Int}) where C = MCI{C}(V, A) # is this really needed?

indices(M::MCI) = collect(1:size(M.A, 2))

function ambient_dim(M::MCI)
    return size(M.A, 1)
end

function rank!(M::MCI, inds::Vector{Int})
    return get(M.rank_cache, inds) do
        rki = rk(M.V[:, inds])
        M.rank_cache[inds] = rki
        rki
    end
end

# --- CellTable --- #

struct CellTable
    cell::MixedCell
    walls::Vector{Hyperplane}
    i_min::Int
    t_min::OrderDual
end

function CellTable(m::MixedCell,
                   M::MCI,
                   p0::DualVector{Int},
                   p1::DualVector{Int},
                   t_curr::OrderDual=zero(OrderDual))

    walls, i_min, t_min = compute_active_walls!(m, M, p0, p1, t_curr)
    return CellTable(m, walls, i_min, t_min)
end

Base.:(==)(a::CellTable, b::CellTable) = a.cell == b.cell
Base.hash(a::CellTable, h::UInt) = Base.hash(a.cell, h)

Base.isless(a::CellTable, b::CellTable) = Base.isless(a.t_min, b.t_min)

is_finished(mtbl::CellTable) = iszero(mtbl.i_min)

# --- WalkData --- # 

mutable struct WalkData{C}
    M::MCI{C}
    cells::MinExtractor{CellTable}
    finished_cells::Set{MixedCellInds}
    p0::DualVector{Int}
    p1::DualVector{Int}
end

function WalkData(M::MCI,
                  initial_mixed_cells::Vector{MixedCell},
                  p0::DualVector{Int},
                  p1::DualVector{Int})

    cells = MinExtractor{CellTable}()
    for m in initial_mixed_cells
        push!(cells, CellTable(m, M, p0, p1))
    end
    return WalkData(M, cells, Set{MixedCellInds}(), p0, p1)
end

# --- ELimination --- #

mutable struct ElimData
    A::Matrix{Int}
    V::Matrix{FqFieldElem}
    shift::Int
    current_covec::Vector{Int}
    current_lift::Vector{Int}
    current_ms::Vector{MixedCellInds}

    function ElimData(A::Matrix{Int}, V::Matrix{C}) where C
        n = size(V, 1) - 1
        k = size(A, 1) - n - 1
        w = rand(-10:10, k + 1)
        Aw = vcat(A[1:n, :], permutedims(vcat(zeros(Int, n), w)) * A)
        shft = min(minimum(i -> Aw[end, i], 1:size(Aw, 2)) - 1, 0) 
        A_shift = copy(Aw)
        A_shift[end, :] = repeat([shft], 1, size(Aw, 2))
        Vp = reduce_mod_rand_prime(V)
        p, ms = with_logger(NullLogger()) do
            mixed_subdivision(hcat(Aw, A_shift), hcat(Vp, Vp),
                              rand(-LSIZE:LSIZE, 2*size(A, 2)))
        end
        return new(A, Vp, shft, w, [pi.eps for pi in p], ms)
    end
end

mutable struct ElimDataDeform
    A::Matrix{Int}
    Ap_deform::Matrix{Int}
    V::Matrix{FqFieldElem}
    Vp::Matrix{FqFieldElem}
    current_lift::DualVector{Int}
    current_ms::Vector{MixedCellInds}

    function ElimDataDeform(A::Matrix{Int}, V::Matrix{C}) where C
        n = size(V, 1) - 1
        VP = reduce_mod_rand_prime(V)
        FF = parent(first(VP))
        Vp = FF.(rand(-1000:1000, n, n + 1)) * VP
        Ap_deform = 5000*A[1:n, :] + rand(-10:10, n, size(A, 2))
        d, ms = mixed_subdivision(Ap_deform, Vp)
        current_lift = [DualNumber(di.eps, 0) for di in d]

        return new(A, Ap_deform, VP, Vp, current_lift, ms)
    end
end
