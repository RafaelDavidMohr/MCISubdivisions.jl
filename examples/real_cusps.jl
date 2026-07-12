using MCISubdivisions
using Oscar
MCIS = MCISubdivisions
using Logging

R, (x, y, z) = QQ[:x, :y, :z]

d = 4
m = ideal(R, gens(R))
f = sum([rand([-1,1])*mon for i in 0:d for mon in gens(m^i)])
F = [f, x*derivative(f, x), x^2*derivative(derivative(f, x), x)]
A, V = MCIS.get_eci_data(F)
d_init = rand(-50000:50000, size(A, 2))
p, ms = mixed_subdivision(A, V, d_init)

function random_lift(l)
    a = randn(l)
    b = [rationalize(x, tol = 0.1) for x in a]
    mult = lcm((denominator).(b))
    return (numerator).(mult * b)
end

p_new = MCIS.DualVector(random_lift(size(A, 2)))
i = 1
max_rr = 0
while true
    ms = with_logger(NullLogger()) do
        MCIS.deform_subdivision(A, V, ms, p, p_new)
    end
    rr = real_root_count(A, QQ.(V), ms)
    println("sample $i, rr $rr")
    if rr > max_rr
        max_rr = rr
        println("NEW MAXIMUM")
    end
    if max_rr == 24
        println("24 real roots reached!")
        break
    end
    p = p_new
    p_new = MCIS.DualVector(random_lift(size(A, 2)))
    i += 1
end
