using MCISubdivisions
using Oscar
MCIS = MCISubdivisions
using Logging

R, (x, y, z) = QQ[:x, :y, :z]

d = 4
m = ideal(R, gens(R))
f = sum([rand(1:1000)*mon for i in 0:d for mon in gens(m^i)])
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

for i in 1:5000
    println(i)
    p_new = MCIS.DualVector(random_lift(size(A, 2)))
    ms = with_logger(NullLogger()) do
        MCIS.deform_subdivision(A, V, ms, p, p_new)
    end
    p = p_new
    mv = sum([MCIS.vol(m, A) for m in ms])
    if mv != 24
        println("wrong mixed volume $(mv)")
        break
    end
    rr = MCIS.real_root_count(ms, A, V)
    if rr == 24
        println("24 real roots!!")
        break
    end
end

