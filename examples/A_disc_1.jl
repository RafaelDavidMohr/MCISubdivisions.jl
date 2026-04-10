using MCISubdivisions
using Oscar

A = [1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1]
A = A[2:end, :]
F = get_A_disc_equations(A)
