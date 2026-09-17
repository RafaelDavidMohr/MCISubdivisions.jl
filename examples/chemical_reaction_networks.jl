using MCISubdivisions
using Oscar
using VerticalRootCounts
using Distributed

# ensure that you are in the same directory as this file

MCIS = MCISubdivisions

function specialize(F::Vector{<:MPolyRingElem}, choice_of_parameters::Vector{<:Union{Int, RingElem}})
    Kax = parent(first(F))
    Ka = coefficient_ring(Kax)
    K = base_ring(Ka)
    Kx, x = polynomial_ring(K, symbols(Kax))
    phi = hom(Kax, Kx, c -> evaluate(c, choice_of_parameters), x)
    return phi.(F)
end

# generic_root_count can take arbitrarily much longer than MCIS.mixed_volume,
# so we run it on a separate worker process which gets killed as soon as the
# computation exceeds TIMEOUT_FACTOR times the time the corresponding mixed
# volume computation took.

const TIMEOUT_FACTOR = 50

const grc_worker = Ref(0) # pid of the worker process, 0 if none is running
const grc_worker_ospid = Ref(0) # operating system pid of that process

function start_grc_worker!()
    grc_worker[] = only(addprocs(1))
    grc_worker_ospid[] = Int(remotecall_fetch(getpid, grc_worker[]))
    Distributed.remotecall_eval(Main, grc_worker[], quote
        using Oscar
        using VerticalRootCounts
        # compile generic_root_count on a tiny example, so that compilation
        # time is not counted towards the timeouts below
        let (C, M, L, _) = multisite_phosphorylation_matrices(1)
            generic_root_count(AugmentedVerticalSystem(C, M, L))
        end
    end)
    return grc_worker[]
end

function kill_grc_worker!()
    if grc_worker[] != 0
        rmprocs(grc_worker[], waitfor = 0)
        # a worker in the middle of a long computation never gets around to
        # processing the exit message sent by rmprocs, so kill it for good
        run(ignorestatus(`kill -9 $(grc_worker_ospid[])`), wait = false)
    end
    grc_worker[] = 0
    grc_worker_ospid[] = 0
    return nothing
end

"""
    timed_generic_root_count(timeout, system::Expr)

Build the `AugmentedVerticalSystem` described by the expression `system` on a
worker process, run `generic_root_count` on it there and return the tuple
`(time, count)`. If the computation does not finish within `timeout` seconds,
the worker process is killed and `(timeout, nothing)` is returned.
"""
function timed_generic_root_count(timeout::Real, system::Expr)
    grc_worker[] == 0 && start_grc_worker!()
    result = Ref{Any}(nothing)
    # a worker in the middle of a computation answers no questions about its
    # state, so we have to wait for it in a task of our own
    computation = @async try
        result[] = remotecall_fetch(Core.eval, grc_worker[], Main, quote
            let F = $system
                tim = @elapsed grc = generic_root_count(F)
                (tim, grc.count)
            end
        end)
    catch err
        result[] = err
    end
    started = time()
    while !istaskdone(computation) && time() - started < timeout
        sleep(0.1)
    end
    if !istaskdone(computation)
        @info "generic_root_count timed out after $timeout seconds"
        kill_grc_worker!() # a fresh worker is started for the next example
        return (Float64(timeout), nothing)
    end
    result[] isa Tuple && return result[]
    # the worker died, most likely because it ran out of memory
    result[] isa Distributed.ProcessExitedException || throw(result[])
    elapsed = time() - started
    @info "the worker computing generic_root_count died after $elapsed seconds"
    kill_grc_worker!()
    return (elapsed, nothing)
end

"""
    timed_generic_root_count(timeout, A, V)

Time out `generic_root_count` on the vertical system given by the exponent
matrix `A` and the coefficient matrix `V` of an ECI system.
"""
function timed_generic_root_count(timeout::Real, A::Matrix{Int}, V::Matrix)
    # the exponents of the ECI systems may be negative; multiplying the system
    # by a monomial changes neither its roots in the torus nor its mixed
    # volume, so shift the exponents of each variable to be non-negative
    M = A .- minimum(A, dims = 2)
    return timed_generic_root_count(timeout, quote
        AugmentedVerticalSystem(matrix(QQ, $V), matrix(ZZ, $M))
    end)
end

my_timings_and_results = Dict{Int, Tuple{Float64, Int}}()
# a count of `nothing` records a timed out computation, the time is then the
# timeout that was exceeded
vrc_timings_and_results = Dict{Int, Tuple{Float64, Union{Int, Nothing}}}()

# example 1
A, V = MCIS.eci_from_odebase("./odebase1.txt", "./odebase1_constraints.txt")
tim = @elapsed mv = MCIS.mixed_volume(A, V)
my_timings_and_results[1] = (tim, mv)

vrc_timings_and_results[1] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), A, V)

# Gro+16
include("./steady_state.jl")
A, V = MCIS.get_eci_data(target_system)
tim = @elapsed mv = MCIS.mixed_volume(A, V)
my_timings_and_results[2] = (tim, mv)

vrc_timings_and_results[2] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), A, V)

# example 2
A, V = MCIS.eci_from_odebase("./odebase2.txt", "./odebase2_constraints.txt")
tim = @elapsed mv = MCIS.mixed_volume(A, V)
my_timings_and_results[3] = (tim, mv)

vrc_timings_and_results[3] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), A, V)

# example 3
A, V = MCIS.eci_from_odebase("./odebase3.txt", "./odebase3_constraints.txt")
tim = @elapsed mv = MCIS.mixed_volume(A, V)
my_timings_and_results[4] = (tim, mv)

vrc_timings_and_results[4] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), A, V)

# example 4
A, V = MCIS.eci_from_odebase("./odebase4.txt", "./odebase4_constraints.txt")
tim = @elapsed mv = MCIS.mixed_volume(A, V)
my_timings_and_results[5] = (tim, mv)

vrc_timings_and_results[5] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), A, V)

# k-site phosphorylation networks

C, M, L, _ = multisite_phosphorylation_matrices(7);
F = AugmentedVerticalSystem(C, M, L).system;
number_of_parameters = ngens(coefficient_ring(parent(first(F))))
target_parameters = rand(1:100, number_of_parameters);
eqns = specialize(F, target_parameters);
tim = @elapsed mv = MCIS.mixed_volume(eqns)
my_timings_and_results[6] = (tim, mv)

vrc_timings_and_results[6] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), quote
    C, M, L, _ = multisite_phosphorylation_matrices(7)
    AugmentedVerticalSystem(C, M, L)
end)

C, M, L, _ = multisite_phosphorylation_matrices(9);
F = AugmentedVerticalSystem(C, M, L).system;
number_of_parameters = ngens(coefficient_ring(parent(first(F))))
target_parameters = rand(1:100, number_of_parameters);
eqns = specialize(F, target_parameters);
tim = @elapsed mv = MCIS.mixed_volume(eqns)
my_timings_and_results[7] = (tim, mv)

vrc_timings_and_results[7] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), quote
    C, M, L, _ = multisite_phosphorylation_matrices(9)
    AugmentedVerticalSystem(C, M, L)
end)

C, M, L, _ = multisite_phosphorylation_matrices(11);
F = AugmentedVerticalSystem(C, M, L).system;
number_of_parameters = ngens(coefficient_ring(parent(first(F))))
target_parameters = rand(1:10000, number_of_parameters);
eqns = specialize(F, target_parameters);
tim = @elapsed mv = MCIS.mixed_volume(eqns)
my_timings_and_results[8] = (tim, mv)

vrc_timings_and_results[8] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), quote
    C, M, L, _ = multisite_phosphorylation_matrices(11)
    AugmentedVerticalSystem(C, M, L)
end)

C, M, L, _ = multisite_phosphorylation_matrices(13);
F = AugmentedVerticalSystem(C, M, L).system;
number_of_parameters = ngens(coefficient_ring(parent(first(F))))
target_parameters = rand(1:10000, number_of_parameters);
eqns = specialize(F, target_parameters);
tim = @elapsed mv = MCIS.mixed_volume(eqns)
my_timings_and_results[9] = (tim, mv)

vrc_timings_and_results[9] = timed_generic_root_count(max(TIMEOUT_FACTOR * tim, 60 * 5), quote
    C, M, L, _ = multisite_phosphorylation_matrices(13)
    AugmentedVerticalSystem(C, M, L)
end)

kill_grc_worker!()
