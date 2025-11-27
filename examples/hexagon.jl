using MCISubdivisions
using Oscar

R, (x, y) = QQ[:x, :y]
r = 1 + y
g = x + x^2*y
b = x^2*y^2 + x*y^2
F = [r + g + b, r + 2*g + 3*b]
