using MCISubdivisions
using Oscar

MCIS = MCISubdivisions

# ensure that you are in the correct directory

# example 1
A, V = MCIS.eci_from_odebase("./odebase1.txt", "./odebase1_constraints.txt")

# example 2
A, V = MCIS.eci_from_odebase("./odebase2.txt", "./odebase2_constraints.txt")

# example 3
A, V = MCIS.eci_from_odebase("./odebase3.txt", "./odebase3_constraints.txt")

# example 4
A, V = MCIS.eci_from_odebase("./odebase4.txt", "./odebase4_constraints.txt")
