using MCISubdivisions
using Oscar

(a10, a11, a12, a13, a14, a20, a21, a22, a23, a24, a30, a31, a32, a33, a34, a40, a41, a42, a43, a44) = rand(1:100,20)
R, (u1, u2, v1, v2, z11, z12, z21, z22, z31, z32, z41, z42, y1, y2) = QQ[:u1, :u2, :v1, :v2, :z11, :z12, :z21, :z22, :z31, :z32, :z41, :z42, :y1, :y2]

u_variables = [u1, u2]
v_variables = [v1, v2]

linear_part = [a10 + a11*u1 + a12*v1 + a13*z11 + a14*z12,
               a20 + a21*u1 + a22*v1 + a23*z21 + a24*z22,
               a30 + a31*u2 + a32*v2 + a33*z31 + a34*z32,
               a40 + a41*u2 + a42*v2 + a43*z41 + a44*z42]

z_variables = [z11, z12, z21, z22, z31, z32, z41, z42]
z_expressions = [u1*y1, u1*y2, v1*y1, v1*y2, u2*y1, u2*y2, v2*y1, v2*y2]

y_variables = [y1, y2]
y_expressions = [u1^2 + v1^2, u2^2 + v2^2]

F = vcat(linear_part, z_variables - z_expressions, y_variables-y_expressions)
