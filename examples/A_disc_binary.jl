using MCISubdivisions
using Oscar

k = 5
A = permutedims(collect(0:k))
F = get_A_disc_equations(A)
