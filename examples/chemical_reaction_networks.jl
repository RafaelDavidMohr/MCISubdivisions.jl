using MCISubdivisions
using Oscar
using VerticalRootCounts

MCIS = MCISubdivisions

function specialize(F::Vector{<:MPolyRingElem}, choice_of_parameters::Vector{<:Union{Int, RingElem}})
    Kax = parent(first(F))
    Ka = coefficient_ring(Kax)
    K = base_ring(Ka)
    Kx, x = polynomial_ring(K, symbols(Kax))
    phi = hom(Kax, Kx, c -> evaluate(c, choice_of_parameters), x)
    return phi.(F)
end

# ensure that you are in the correct directory

# example 1
A, V = MCIS.eci_from_odebase("./odebase1.txt", "./odebase1_constraints.txt")
MCIS.mixed_volume(A, V)

# Gro+16
include("./steady_state.jl")
A, V = MCIS.get_eci_data(target_system)
MCIS.mixed_volume(A, V)

# example 2
A, V = MCIS.eci_from_odebase("./odebase2.txt", "./odebase2_constraints.txt")
MCIS.mixed_volume(A, V)

# example 3
A, V = MCIS.eci_from_odebase("./odebase3.txt", "./odebase3_constraints.txt")
MCIS.mixed_volume(A, V)

# example 4
A, V = MCIS.eci_from_odebase("./odebase4.txt", "./odebase4_constraints.txt")
MCIS.mixed_volume(A, V)

# k-site phosphorylation networks
C, M, L, _ = multisite_phosphorylation_matrices(9);
F = AugmentedVerticalSystem(C, M, L).system;
number_of_parameters = ngens(coefficient_ring(parent(first(F))))
target_parameters = rand(1:100, number_of_parameters);
eqns = specialize(F, target_parameters);
MCIS.mixed_volume(eqns)

C, M, L, _ = multisite_phosphorylation_matrices(11);
F = AugmentedVerticalSystem(C, M, L).system;
number_of_parameters = ngens(coefficient_ring(parent(first(F))))
target_parameters = rand(1:10000, number_of_parameters);
eqns = specialize(F, target_parameters);
MCIS.mixed_volume(eqns)

C, M, L, _ = multisite_phosphorylation_matrices(13);
F = AugmentedVerticalSystem(C, M, L).system;
number_of_parameters = ngens(coefficient_ring(parent(first(F))))
target_parameters = rand(1:10000, number_of_parameters);
eqns = specialize(F, target_parameters);
MCIS.mixed_volume(eqns)
