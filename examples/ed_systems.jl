using Oscar
using MCISubdivisions

MCIS = MCISubdivisions

function ed_system(eqns)
    R = parent(first(eqns))
    n = ngens(R)
    FF = base_ring(R)
    S, z, x, u = polynomial_ring(FF, "z" => 1:length(eqns), "x" => 1:n, "u" => 1:n)
    phi = hom(R, S, x)
    F = phi.(eqns)
    return vcat(F, (x .* (u .- x)) - sum([[zz * xi * derivative(phi(f), xi) for xi in x] for (zz, f) in zip(z, eqns)]))
end

R, x = polynomial_ring(GF(65521), ["x$i" for i in 1:3])

supp = [rand(0:4, 3) for _ in 1:9];
push!(supp, zeros(Int, 3));
F = [sum([rand(1:65520) * prod(x .^ Int.(v)) for v in supp])
     for _ in 1:2];
eqns = ed_system(F);

supp = [rand(0:4, 3) for _ in 1:13];
push!(supp, zeros(Int, 3));
F = [sum([rand(1:65520) * prod(x .^ Int.(v)) for v in supp])
     for _ in 1:2];
eqns = ed_system(F);

supp = [rand(0:4, 3) for _ in 1:17];
push!(supp, zeros(Int, 3));
F = [sum([rand(1:65520) * prod(x .^ Int.(v)) for v in supp])
     for _ in 1:2];
eqns = ed_system(F);

supp = [rand(0:4, 3) for _ in 1:21];
push!(supp, zeros(Int, 3));
F = [sum([rand(1:65520) * prod(x .^ Int.(v)) for v in supp])
     for _ in 1:2];
eqns = ed_system(F);

R, x = polynomial_ring(GF(65521), ["x$i" for i in 1:4])

supp = [rand(0:4, 4) for _ in 1:5];
push!(supp, zeros(Int, 4));
F = [sum([rand(1:65520) * prod(x .^ Int.(v)) for v in supp])
     for _ in 1:3];
eqns = ed_system(F);

supp = [rand(0:4, 4) for _ in 1:7];
push!(supp, zeros(Int, 4));
F = [sum([rand(1:65520) * prod(x .^ Int.(v)) for v in supp])
     for _ in 1:3];
eqns = ed_system(F);

supp = [rand(0:4, 4) for _ in 1:9];
push!(supp, zeros(Int, 4));
F = [sum([rand(1:65520) * prod(x .^ Int.(v)) for v in supp])
     for _ in 1:3];
eqns = ed_system(F);

# evaluate with examples above
A, V = MCIS.get_eci_data(eqns);
@elapsed E = MCIS.ElimData(A, V) # initial subdivision
covec = rand(-100:100, size(A, 1) - size(V, 1) + 1);
@elapsed MCIS.elim_vertex!(E, covec)
