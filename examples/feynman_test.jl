using MCISubdivisions
using Oscar

MCIS = MCISubdivisions

function dehomogenize(eqns, i)
    R = parent(first(eqns))
    S, x = polynomial_ring(base_ring(R), symbols(R)[1:end .!= i])
    phi = hom(R, S, vcat(x[1:i-1], [one(S)], x[i:end]))
    return phi.(eqns)
end

# debox
R, (a_1,a_2,a_3,a_4,a_5,M,m,s,t) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:5], ["M", "m", "s", "t"])) 
U_debox = a_1*a_4 + a_2*a_4 + a_3*a_4 + a_1*a_5 + a_2*a_5 + a_3*a_5 + a_4*a_5;
F0 = M*a_1*a_2*a_4 + s*a_1*a_3*a_4 + M*a_2*a_3*a_4 + M*a_1*a_2*a_5 + s*a_1*a_3*a_5 + M*a_2*a_3*a_5 + M*a_1*a_4*a_5 + t*a_2*a_4*a_5 + M*a_3*a_4*a_5;
F = F0 - m*U_debox*(a_1+a_2+a_3+a_4+a_5);
a = gens(R)[1:5]
eqns = dehomogenize([ai*derivative(F, ai) for ai in a], 5)

# npltrb
R, (a_1,a_2,a_3,a_4,a_5,a_6,M,m,s) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:6], ["M", "m", "s"])) 
U_npltrb = a_1*a_2 + a_1*a_3 + a_1*a_4 + a_2*a_4 + a_3*a_4 + a_2*a_5 + a_3*a_5 + a_4*a_5 + a_1*a_6 + a_2*a_6 + a_3*a_6 + a_5*a_6;
F0 = M*a_1*a_2*a_3 + M*a_1*a_2*a_4 + s*a_1*a_3*a_4 + M*a_2*a_3*a_4 + s*a_1*a_2*a_5 + s*a_1*a_3*a_5 + M*a_2*a_3*a_5 + s*a_1*a_4*a_5 + M*a_2*a_4*a_5 + M*a_1*a_3*a_6 + M*a_2*a_3*a_6 + M*a_1*a_4*a_6 + M*a_2*a_4*a_6 + M*a_3*a_4*a_6 + s*a_1*a_5*a_6 + s*a_2*a_5*a_6 + M*a_3*a_5*a_6 + M*a_4*a_5*a_6;
F = F0 - m*U_npltrb*(a_1+a_2+a_3+a_4+a_5+a_6);
a = gens(R)[1:6]
eqns = dehomogenize([ai*derivative(F, ai) for ai in a], 6)
