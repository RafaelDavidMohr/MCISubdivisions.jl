using MCISubdivisions
using Oscar

MCIS = MCISubdivisions

function dehomogenize(eqns, i)
    R = parent(first(eqns))
    S, x = polynomial_ring(base_ring(R), symbols(R)[1:end .!= i])
    phi = hom(R, S, vcat(x[1:i-1], [one(S)], x[i:end]))
    return phi.(eqns)
end

# npltrb
R, (a_1,a_2,a_3,a_4,a_5,a_6,M,m,s) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:6], ["M", "m", "s"])) 
U_npltrb = a_1*a_2 + a_1*a_3 + a_1*a_4 + a_2*a_4 + a_3*a_4 + a_2*a_5 + a_3*a_5 + a_4*a_5 + a_1*a_6 + a_2*a_6 + a_3*a_6 + a_5*a_6;
F0 = M*a_1*a_2*a_3 + M*a_1*a_2*a_4 + s*a_1*a_3*a_4 + M*a_2*a_3*a_4 + s*a_1*a_2*a_5 + s*a_1*a_3*a_5 + M*a_2*a_3*a_5 + s*a_1*a_4*a_5 + M*a_2*a_4*a_5 + M*a_1*a_3*a_6 + M*a_2*a_3*a_6 + M*a_1*a_4*a_6 + M*a_2*a_4*a_6 + M*a_3*a_4*a_6 + s*a_1*a_5*a_6 + s*a_2*a_5*a_6 + M*a_3*a_5*a_6 + M*a_4*a_5*a_6;
F = F0 - m*U_npltrb*(a_1+a_2+a_3+a_4+a_5+a_6);
a = gens(R)[1:6]
eqns = dehomogenize([ai*derivative(F, ai) for ai in a], 6)

# pltrb
R, (a_1,a_2,a_3,a_4,a_5,a_6,M,m,s) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:6], ["M", "m", "s"])) 
U_pltrb = a_1*a_2 + a_1*a_3 + a_1*a_4 + a_2*a_5 + a_3*a_5 + a_4*a_5 + a_1*a_6 + a_2*a_6 + a_3*a_6 + a_4*a_6 + a_5*a_6;
F0 = M*a_1*a_2*a_3 + s*a_1*a_2*a_4 + M*a_1*a_3*a_4 + s*a_1*a_2*a_5 + s*a_1*a_3*a_5 + M*a_2*a_3*a_5 + s*a_1*a_4*a_5 + s*a_2*a_4*a_5 + M*a_3*a_4*a_5 + M*a_1*a_3*a_6 + M*a_2*a_3*a_6 + s*a_1*a_4*a_6 + s*a_2*a_4*a_6 + M*a_3*a_4*a_6 + s*a_1*a_5*a_6 + s*a_2*a_5*a_6 + M*a_3*a_5*a_6;
F = F0 - m*U_pltrb*(a_1+a_2+a_3+a_4+a_5+a_6);
a = gens(R)[1:6]
eqns = dehomogenize([ai*derivative(F, ai) for ai in a], 6)

# dbox
R, (a_1,a_2,a_3,a_4,a_5,a_6,a_7,M,m,s,t) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:7], ["M", "m", "s", "t"])) 
U_dbox = a_1*a_3 + a_2*a_3 + a_1*a_4 + a_2*a_4 + a_1*a_5 + a_2*a_5 + a_3*a_6 + a_4*a_6 + a_5*a_6 + a_1*a_7 + a_2*a_7 + a_3*a_7 + a_4*a_7 + a_5*a_7 + a_6*a_7;
F0 = M*a_1*a_2*a_3 + M*a_1*a_2*a_4 + M*a_1*a_3*a_4 + M*a_2*a_3*a_4 + M*a_1*a_2*a_5 + s*a_1*a_3*a_5 + s*a_2*a_3*a_5 + M*a_1*a_4*a_5 + M*a_2*a_4*a_5 + M*a_1*a_3*a_6 + s*a_2*a_3*a_6 + M*a_1*a_4*a_6 + s*a_2*a_4*a_6 + M*a_3*a_4*a_6 + M*a_1*a_5*a_6 + s*a_2*a_5*a_6 + s*a_3*a_5*a_6 + M*a_4*a_5*a_6 + M*a_1*a_2*a_7 + M*a_1*a_3*a_7 + t*a_1*a_4*a_7 + M*a_2*a_4*a_7 + M*a_3*a_4*a_7 + M*a_1*a_5*a_7 + s*a_2*a_5*a_7 + s*a_3*a_5*a_7 + M*a_4*a_5*a_7 + M*a_1*a_6*a_7 + s*a_2*a_6*a_7 + s*a_3*a_6*a_7 + M*a_4*a_6*a_7;
F = F0 - m*U_dbox*(a_1+a_2+a_3+a_4+a_5+a_6+a_7);
a = gens(R)[1:7]
eqns = dehomogenize([ai*derivative(F, ai) for ai in a], 7)

# tdetri
R, (a_1,a_2,a_3,a_4,a_5,M,m,s) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:5], ["M", "m", "s"])) 
U_tdetri = a_1*a_2*a_3 + a_1*a_2*a_4 + a_1*a_3*a_4 + a_2*a_3*a_4 + a_1*a_3*a_5 + a_2*a_3*a_5 + a_1*a_4*a_5 + a_2*a_4*a_5;
F0 = s*a_1*a_2*a_3*a_4 + M*a_1*a_2*a_3*a_5 + M*a_1*a_2*a_4*a_5 + M*a_1*a_3*a_4*a_5 + M*a_2*a_3*a_4*a_5;
F = F0 - m*U_tdetri*(a_1+a_2+a_3+a_4+a_5);
a = gens(R)[1:5]
eqns = dehomogenize([ai*derivative(F, ai) for ai in a], 5)
eqns = dehomogenize(eqns, 3)

# debox
R, (a_1,a_2,a_3,a_4,a_5,M,m,s,t) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:5], ["M", "m", "s", "t"])) 
U_debox = a_1*a_4 + a_2*a_4 + a_3*a_4 + a_1*a_5 + a_2*a_5 + a_3*a_5 + a_4*a_5;
F0 = M*a_1*a_2*a_4 + s*a_1*a_3*a_4 + M*a_2*a_3*a_4 + M*a_1*a_2*a_5 + s*a_1*a_3*a_5 + M*a_2*a_3*a_5 + M*a_1*a_4*a_5 + t*a_2*a_4*a_5 + M*a_3*a_4*a_5;
F = F0 - m*U_debox*(a_1+a_2+a_3+a_4+a_5);
a = gens(R)[1:5]
eqns = [ai*derivative(F, ai) for ai in a];
eqns = dehomogenize(eqns, 5)

# acn
R, (a_1,a_2,a_3,a_4,a_5,M,m,s,t) = polynomial_ring(QQ, vcat(["a_$i" for i in 1:5], ["M", "m", "s", "t"])) 
U_acn = a_5*(a_1+a_2+a_3+a_4) + a_1*(a_2+a_3)+a_4*(a_2+a_3)
F1 = a_1*a_2*a_5
F2 = a_2*a_3*(a_1+a_4+a_5)
F3 = a_3*a_4*a_5
F4 = a_1*a_4*(a_2+a_3+a_5)
F12 = a_1*a_3*a_5
F23 = a_2*a_4*a_5
F13 = 0 
F = F12*s + F23*t + M*(F1+F2+F3+F4) - U_acn*m*(a_1+a_2+a_3+a_4+a_5)
a = gens(R)[1:5]
# eqns = [ai*derivative(F, ai) for ai in a];
eqns = [ai*derivative(F, ai) for ai in a];
eqns = dehomogenize(eqns, 1)
eqns = dehomogenize(eqns, 8)


