using Oscar
using MCISubdivisions

MCIS = MCISubdivisions

function critical_value_system(supp, lx, ly)
    _, y, z, x = polynomial_ring(GF(65521), "y" => 1:ly, "z" => 1:(ly-1), "x" => 1:lx)
    xy = vcat(x, y)
    F = [sum([rand(1:65520) * prod(xy .^ s) for s in supp])]
    for i in 1:(ly-1)
        push!(F, z[i]*sum([rand(1:65520) * prod(xy .^ s) for s in supp]))
    end
    for i in 1:ly
        push!(F, sum([y[i] * derivative(f, y[i]) for f in F]))
    end
    return F
end

function test_polytope(rhs)
    n = length(rhs)
    A = [j == i ? 1 : 0 for j in 1:n, i in 1:n]
    A = vcat(A, -A)
    A = vcat(A, ones(Int, 1, n))
    return polyhedron(A, vcat(rhs, zeros(Int, n), [Int(ceil(sum(rhs) / 2))]))
end

function ed_thing(eqns)
    R = parent(first(eqns))
    n = ngens(R)
    FF = base_ring(R)
    S, z, x, u = polynomial_ring(FF, "z" => 1:length(eqns), "x" => 1:n, "u" => 1:n)
    phi = hom(R, S, x)
    F = phi.(eqns)
    return vcat(F, (x .* (u .- x)) - z .* x .* [derivative(phi(f), xi) for xi in x])
end
    

# use with three lines of code below
P = test_polytope([5,4,3,3]);
supp = [Int.(s) for s in vertices(P)];
F = critical_value_system(supp, 2, 2);
Q = get_eliminant_polytope(F);

P = test_polytope([6,4,3,3,3])
supp = [Int.(s) for s in vertices(P)];
F = critical_value_system(supp, 2, 3);
Q = get_eliminant_polytope(F);



