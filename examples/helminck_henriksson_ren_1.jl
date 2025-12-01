using MCISubdivisions
using Oscar

R, (x1, x2) = QQ[:x1, :x2]
t = rand(1:100, 5) .* [x1^2, x2^2, x1, x2, R(1)]
F = [t[1] + t[2] + t[3] + t[4] + t[5],
     3*t[1] + 3*t[2] + 5*t[3] + 7*t[4] + 11*t[5]]
