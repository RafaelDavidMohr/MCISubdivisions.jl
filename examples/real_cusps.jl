using MCISubdivisions
using Oscar
MCIS = MCISubdivisions
using Logging

R, (x, y, z) = QQ[:x, :y, :z]

d = 4
m = ideal(R, gens(R))
f = sum([mon for i in 0:d for mon in gens(m^i)])
F = [f, derivative(f, x), derivative(derivative(f, x), x)]
A, V = MCIS.get_eci_data(F)
rand_mix = matrix(QQ, (QQ).(rand(-1000:1000, size(V, 1), size(V, 1))))
V = Matrix(rand_mix * matrix(QQ, V))
p, ms = mixed_subdivision(A, V)
sum([MCIS.vol(m, A) for m in ms])

function random_lift(l)
    a = randn(l)
    b = [rationalize(x, tol = 0.01) for x in a]
    mult = lcm((denominator).(b))
    return (numerator).(mult * b)
end

t_cross = MCIS.dual_zero(GF(2))
p_new = MCIS.DualVector(random_lift(size(A, 2)))
i = 1
while true
    println("sample $i")
    reached, t_cross, ms = with_logger(NullLogger()) do
        MCIS.deform_subdivision(A, V, ms, p, p_new, 24)
    end
    if reached
        println("24 real roots!!")
        break
    end
    p = p_new
    p_new = MCIS.DualVector(random_lift(size(A, 2)))
    i += 1
end
