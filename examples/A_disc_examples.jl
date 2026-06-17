using MCISubdivisions
using Oscar

# Example 1
A = [1 1 1 1 0 0 0 0; 0 0 0 0 1 1 1 1; 2 3 5 7 11 13 17 19; 19 17 13 11 7 5 3 2]
A = A[2:end, :]
F = get_A_disc_equations(A)

# Example 2
A = [1 1 1 1 1 0 0 0 0 0; 0 0 0 0 0 1 1 1 1 1; 2 3 5 7 11 13 17 19 21 23; 23 21 19 17 13 11 7 5 3 2]
A = A[2:end, :]
F = get_A_disc_equations(A)

# Example 3
A = [1 0 1 1 1 0 1 0 1 0 1 1; 1 1 1 1 0 0 0 1 0 1 1 0; 1 0 1 0 0 0 1 0 0 1 0 1; 1 0 0 0 1 0 1 1 0 1 1 0]
F = get_A_disc_equations(A)

# Δ_1 × Δ_1 × Δ_1
vs = vertices(simplex(1) * simplex(1) * simplex(1))
A = Int.(hcat(vs...))
F = get_A_disc_equations(A)

# Δ_1 × Δ_1 × Δ_2
vs = vertices(simplex(1) * simplex(1) * simplex(2))
A = Int.(hcat(vs...))
F = get_A_disc_equations(A)
