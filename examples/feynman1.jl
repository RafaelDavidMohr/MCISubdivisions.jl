using MCISubdivisions
using Oscar

MCIS = MCISubdivisions

R, (a1,a2,a3,a4,u,m1,m2,m3,m4,M1,M2,M3,M4,s,t) = polynomial_ring(QQ, vcat(["a$i" for i in 1:4], ["u"], ["m$i" for i in 1:4], ["M$i" for i in 1:4], ["s","t"]))

U = a1 + a2 + a3 + a4 
F1 = a1*a2
F2 = a2*a3
F3 = a3*a4
F4 = a1*a4
F12 = a1*a3
F23 = a2*a4

F = (s*F12 + t*F23 + M1*F1 + M2*F2 + M3*F3 + M4*F4 - U*(m1*a1+m2*a2+m3*a3+m4*a4));

a = [a1,a2,a3,a4]
eqns = vcat([ai*derivative(F, ai) for ai in a], [1 - u*prod(a)*U])

function dehomogenize(eqns, i)
    R = parent(first(eqns))
    S, x = polynomial_ring(base_ring(R), symbols(R)[1:end .!= i])
    phi = hom(R, S, vcat(x[1:i-1], [one(S)], x[i:end]))
    return phi.(eqns)
end

eqns_dehom = dehomogenize(eqns, 4)
