# --- mixed cells --- #

function MixedCell(inds::MixedCellInds, M::MCI,
                   known::Dict{Int, Vector{Int}}=Dict{Int, Vector{Int}}())

    for ind in inds
        sort!(ind)
    end
    return MixedCell(inds, find_loc_indices!(M, inds, known))
end

function find_loc_indices!(M::MCI, inds::MixedCellInds,
                           known::Dict{Int, Vector{Int}})

    loc_inds = Vector{Int}[]
    rem_inds = indices(M)
    all_m_inds = Int[]
    rk = 0
    k = length(inds)
    @inbounds for i in 1:length(inds)
        rem_inds = setdiff(rem_inds, inds[i])
        all_m_inds = vcat(all_m_inds, inds[i])
        rk += length(inds[i]) - 1
        loc_inds_i = if i != k
            get(known, i) do
                find_loc_indices!(M, all_m_inds, rem_inds, rk)
            end
        else
            copy(rem_inds)
        end
        push!(loc_inds, loc_inds_i)
        setdiff!(rem_inds, loc_inds_i)
    end
    return loc_inds
end

function find_loc_indices!(M::MCI, inds::Vector{Int}, test_inds::Vector{Int},
                           V_rank::Int) 

    res = Int[]
    @inbounds for i in test_inds
        if rank!(M, vcat(inds, [i])) == V_rank
            push!(res, i)
        end
    end
    return sort(res)
end

function vol(m::MixedCell, A::Matrix{Int})
    return vol(m.inds, A)
end

function vol(m::MixedCellInds, A::Matrix{Int})
    mat = hcat([linear_span(A, S) for S in m]...)
    lu = LinearAlgebra.lu(mat)
    return round(Int, abs(LinearAlgebra.det(lu)))
end

# naive computation using msolve
function real_root_count(m::MixedCellInds, A::Matrix{Int}, V::Matrix{QQFieldElem})
    R, x = polynomial_ring(QQ, ["x$i" for i in 1:size(A,1)])
    monA(i) = prod(x .^ A[:, i])
    mmat = Matrix(echelon_form(matrix(QQ, V)[:, vcat(m...)]))
    eqns = QQMPolyRingElem[]
    for inds in m
        for l in 1:length(inds)-1
            push!(eqns, sum([mmat[l, i] * monA(j) for (i, j) in enumerate(inds)]))
        end
    end
    I = ideal(R, eqns)
    I = saturation(I, ideal(R, [prod(gens(R))]))
    return length(Oscar.real_solutions(I)[1])
end

function cayley_indices(m::MixedCell, j::Int)
    if j == length(m)
        return @inbounds m.loc_inds[j]
    end
    res = @inbounds vcat(m.loc_inds[j], [first(m.inds[j+1])])
    return sort(res)
end

# --- matrix --- #

# returns only first kernel element, hopefully reliable
function krnel_nz_entries(V::Matrix{Int}, indmax::Int;
                          tol = 1e-8)
    k = nullspace(V)[:, 1]
    return findall(j -> abs(k[j]) > tol, 1:indmax)
end

function krnel_nz_entries(V::Matrix{FqFieldElem}, indmax::Int)
    F = parent(first(V))
    K = kernel(matrix(F, V), side = :right)[:, 1]
    return findall(j -> !iszero(K[j]), 1:indmax)
end

rk(V::Matrix{Int}) = rank(V)

function rk(V::Matrix{FqFieldElem})
    F = parent(first(V))
    return rank(matrix(F, V))
end

function reduce_mod_rand_prime(V::Matrix{QQFieldElem})
    p = Hecke.rand_bits_prime(ZZ, 31)
    F = GF(p)
    return [F(numerator(x)) * F(denominator(x))^(-1) for x in V]
end

function reduce_mod_rand_prime(V::Matrix{Int})
    p = Hecke.rand_bits_prime(ZZ, 31)
    F = GF(p)
    return [F(x) for x in V]
end

function linear_span(A::Matrix{C}, inds::Vector{Int}) where C
    a0 = A[:, first(inds)]
    L = Matrix{C}(undef, size(A, 1), length(inds) - 1)
    for (j, i) in enumerate(inds[2:end])
        L[:, j] = A[:, i] - a0
    end
    return L
end

function reduce_to_linear_span(A::Matrix{Int})
    L = linear_span(A, collect(1:size(A, 2)))
    Lt = hermite_form(matrix(ZZ, L), trim = :true)
    return hcat(zeros(Int, size(Lt, 1)), (Int).(Matrix(Lt)))
end

function is_affine_independent(A::Matrix{Int}, inds::MixedCellInds)
    mat = hcat([linear_span(A, S) for S in inds]...)
    return !iszero(round(det(mat))) 
end

function normal_space(A::Matrix{C}, m::MixedCellInds) where C

    n = size(A, 1)
    eqns = Matrix{C}(undef, n, 0)
    for S in m
        a = A[:, first(S)]
        for i in S[2:end]
            eqns = hcat(eqns, A[:, i] - a)
        end
    end

    if C <: Float64
        return nullspace(eqns)
    end

    F = if C <: Int
        QQ
    else
        parent(first(A))
    end

    return permutedims(Matrix(kernel(matrix(F, eqns))))
end

function select_max_weight_columns(A::Matrix, w::Vector)
    column_dict = Dict()
    
    @inbounds for i in 1:size(A, 2)
        col_key = Tuple(A[:, i])
        
        if haskey(column_dict, col_key)
            best_idx, best_weight = column_dict[col_key]
            if w[i] > best_weight
                column_dict[col_key] = (i, w[i])
            end
        else
            column_dict[col_key] = (i, w[i])
        end
    end
    
    return sort([idx for (idx, _) in values(column_dict)])
end

# --- circuits --- #

function partial_sum(inds::Vector{Int}, c::SparseVector{Int, Int})
    return @inbounds sum(c[inds])
end

function partial_sum_sign(inds::Vector{Int}, c::Hyperplane)
    res = partial_sum(inds, c.cfs)
    return res != 0 && (c.sgn ? res > 0 : res < 0) # TODO: check if 0 is allowed here
end

# --- homotopy paths --- #
    
function crossing_val(l0::DualVector, l1::DualVector, c::SparseVector{Int, Int})
    l0d = dot(l0, c)
    l1d = dot(l1, c)
    return l0d, l1d, cross_val_from_dp(l0d, l1d)
end

function cross_val_from_dp(l0d::DualNumber{T}, l1d::DualNumber{T}) where T
    denom = l0d - l1d

    if iszero(l0d.r) && iszero(l1d.r)
        return DualNumber{T}(l0d.eps * denom.eps^(-1), zero(T))
    elseif iszero(denom.r)
        return DualNumber{T}(T(2), zero(T))
    else
        return l0d * inv(denom)
    end
end

# --- other helpers --- #

# find c, d such that c*a - d*b = 0
function zero_coefficients(a::Integer, b::Integer)
    l = lcm(a, b)
    return div(l, a), div(l, b)
end

function zero_coefficients(a::DualNumber{T}, b::DualNumber{T}) where T
    isinvertible(a) && return b*inv(a), one(a)
    isinvertible(b) && return one(a), a*inv(b)
    return DualNumber(b.eps, 0), DualNumber(a.eps, 0)
end

function forgetful_lift(A_size::Int, forget_inds::Vector{Int}, ::T) where T
    d = Vector{T}(undef, A_size)
    @inbounds for i in 1:A_size
        if i in forget_inds
            d[i] = 0
        else
            d[i] = 1
        end
    end
    return DualVector(d)
end

function forgetful_lift(A_size::Int, forget_inds::Vector{Int},
                        d_eps::Vector{T}) where T

    d = Vector{T}(undef, A_size)
    d_eps_new = Vector{T}(undef, A_size)
    cnt = 1
    @inbounds for i in 1:A_size
        if i in forget_inds
            d[i] = 0
            d_eps_new[i] = 0
        else
            d[i] = 1
            d_eps_new[i] = d_eps[cnt]
            cnt += 1
        end
    end
    return [DualNumber(dd, dd_eps) for (dd, dd_eps) in zip(d, d_eps_new)]
end

function get_support(F::Vector{<:MPolyRingElem})
    exps = unique(vcat([collect(exponents(f)) for f in F]...))
    return hcat(exps...)
end

function get_eci_data(F::Vector{<:MPolyRingElem})
    A = get_support(F)
    return A, [coeff(f, A[:, i]) for f in F, i in 1:size(A, 2)]
end

function get_eci_data(F::Vector{QQMPolyRingElem})
    A = get_support(F)
    Vcfs = [coeff(f, A[:, i]) for f in F, i in 1:size(A, 2)]
    V = try
        Int.(Vcfs)
    catch
        "Only integer or finite field coefficients supported."
    end
    return A, V
end

function rand_arr_ff(F::FqField, dims...)
    return (F).(rand(0:characteristic(F)-1, dims...))
end

function gather_mixed_cells(wd::WalkData, max_ind=0::Int)
    iszero(max_ind) && return collect(wd.finished_cells)
    result = MixedCellInds[]
    for minds in wd.finished_cells
        any(S -> any(i -> i > max_ind, S), minds) && continue
        push!(result, minds)
    end
    return result
end

function gather_mixed_cells(wd::WalkData, excluded_inds::Vector{Int})

    isempty(excluded_inds) && return collect(wd.finished_cells)
    ind_map = Dict{Int, Int}()
    rem_inds = sort(setdiff(1:size(wd.M.A, 2), excluded_inds))
    for (j, i) in enumerate(rem_inds)
        ind_map[i] = j
    end
    result = MixedCellInds[]
    for minds in wd.finished_cells
        any(S -> any(i -> i in excluded_inds, S), minds) && continue
        push!(result, [(idx -> ind_map[idx]).(S) for S in minds])
    end

    return result
end

function id_matrix(n)
    return [i == j ? 1 : 0 for i in 1:n, j in 1:n]
end

function get_A_disc_equations(A::Matrix{Int}, dehom=1)

    d, n = size(A)
    A_lift = vcat(A, id_matrix(n))
    A_lift = A_lift[1:size(A_lift, 1) .!= (dehom + d), :]
    R, vars = polynomial_ring(QQ, vcat(["x$i" for i in 1:d], ["z$i" for i in 1:(n-1)]))
    x = vars[1:d]
    z = vars[d+1:end]
    s = sum(prod(map((i,j) -> i^j, vars, A_lift[:,k])) for k in 1:n)
    fs = [s; [v*derivative(s, v) for v in x]]
    return fs
end

function specialize(F::Vector{<:MPolyRingElem},
                    choice_of_parameters::Vector{<:Union{Int, RingElem}})
    
    Kax = parent(first(F))
    Ka = coefficient_ring(Kax)
    K = base_ring(Ka)
    Kx, x = polynomial_ring(K, symbols(Kax))
    phi = hom(Kax, Kx, c -> evaluate(c, choice_of_parameters), x)
    return phi.(F)
end

function eci_from_odebase(ode_filename, constr_filename)
    vrs, prms, rhs = parse_ode_system_manual(ode_filename, constr_filename)
    ss = create_steady_state(vrs, prms, rhs)
    return make_eci(ss)
end

# Functions to parse sage files from odebase.org
function parse_ode_system_manual(ode_filename, constr_filename)
    ode_content = read(ode_filename, String)
    
    ode_rhs_list = String[]
    
    ode_pattern = r"Derivative\([^,]+,\s*t\)\s*,\s*(.+?)\s*\)(?:,| \])"
    
    for m in eachmatch(ode_pattern, ode_content)
        rhs = strip(m.captures[1])
        rhs = replace(rhs, r"\s*\)\s*$" => "")
        push!(ode_rhs_list, "R(" * rhs * ")")
    end
    
    constr_content = read(constr_filename, String)
    
    constr_rhs_list = String[]
    c_params = String[]
    
    constr_pattern = r"Eq\(\s*(.+?)\s*,\s*(.+?)\s*\)"
    
    i = 1
    for m in eachmatch(constr_pattern, constr_content)
        lhs = strip(m.captures[1])
        c_param = "c_$i"
        push!(c_params, c_param)
        push!(constr_rhs_list, "R(" * "($lhs) - $c_param" * ")")
        i += 1
    end
    
    rhs_list = vcat(ode_rhs_list, constr_rhs_list)
    
    all_content = ode_content * " " * constr_content
    var_pattern = r"x\d+"
    var_matches = eachmatch(var_pattern, all_content)
    vars = unique([m.match for m in var_matches])
    x_vars = sort(vars, by=x -> parse(Int, match(r"x(\d+)", x).captures[1]))
    
    param_pattern = r"k\d+"
    param_matches = eachmatch(param_pattern, all_content)
    params = unique([m.match for m in param_matches])
    params_sorted = sort(params, by=k -> parse(Int, match(r"k(\d+)", k).captures[1]))
    
    params_sorted = vcat(params_sorted, c_params)
    
    return x_vars, params_sorted, rhs_list
end

function create_steady_state(x_vars, params, rhs)
    param_ring_expr = quote
        P, $(Expr(:tuple, [Symbol(v) for v in params]...)) = polynomial_ring(QQ, $(params))
        KP = fraction_field(P)
    end
    eval(param_ring_expr)

    poly_ring_expr = quote
        R, $(Expr(:tuple, [Symbol(x) for x in x_vars]...)) = polynomial_ring(KP, $(x_vars))
    end
    eval(poly_ring_expr)

    eqns = [eval(Meta.parse(eqn)) for eqn in rhs]
    return eqns
end

function make_eci(steady_state_system)
    A, V = get_eci_data(steady_state_system)
    A_red = reduce_to_linear_span(A)
    R = parent(first(steady_state_system))
    P = base_ring(R)
    param_evals = rand(1:1000, ngens(P))
    V_qq = [f(param_evals...) for f in V]
    V_qq = rand(-100:100, size(A_red, 1), size(V_qq, 1)) * V_qq
    return A_red, reduce_mod_rand_prime(V_qq)
end

function convert_to_dual_number_vector(p::DualVector)
    return [DualNumber(a, b) for (a, b) in zip(p.r, p.eps)]
end
