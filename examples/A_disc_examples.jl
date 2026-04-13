using MCISubdivisions
using Oscar

A = [1 1 1 1 1 1 1 1; 0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1]
A = A[2:end, :]
F = get_A_disc_equations(A)

A = [1 1 1 1 0 0 0 0; 0 0 0 0 1 1 1 1; 2 3 5 7 11 13 17 19; 19 17 13 11 7 5 3 2]
A = A[2:end, :]
F = get_A_disc_equations(A)

# too large
A = [1 1 1 1 1 1 0 0 0 0 0 0; 0 0 0 0 0 0 1 1 1 1 1 1; 2 3 5 7 11 13 17 19 21 23 25 27; 27 25 23 21 19 17 13 11 7 5 3 2]
A = A[2:end, :]
F = get_A_disc_equations(A)
