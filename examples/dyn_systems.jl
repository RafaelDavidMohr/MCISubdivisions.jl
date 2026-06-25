using MCISubdivisions
using Oscar
using MixedFiberPolytope # https://github.com/RafaelDavidMohr/MixedFiberPolytope

MCIS = MCISubdivisions

function eci_from_dynsys(f, g)
    S = parent(f)
    FF = base_ring(S)
    R, (y, y1, p, q, x, x1, x2) = polynomial_ring(FF, ["y", "y'", "p", "q", "x", "x'", "x''"])
    phi = hom(S, R, [x, y])
    ff = phi(f)
    gg = phi(g)
    return [x1 - ff, y1 - gg,
            x*p - x*derivative(ff, x),
            y*q - y*derivative(ff, y),
            x2 - (x1 * p + y1 * q)]
end

function eci_from_dynsys_for_mfp(f, g) # MixedFiberPolytope eliminates the last block of variables
    S = parent(f)
    FF = base_ring(S)
    R, (x, x1, x2, y, y1, p, q) = polynomial_ring(FF, ["x", "x'", "x''", "y", "y'", "p", "q"])
    phi = hom(S, R, [x, y])
    ff = phi(f)
    gg = phi(g)
    return [x1 - ff, y1 - gg,
            x*p - x*derivative(ff, x),
            y*q - y*derivative(ff, y),
            x2 - (x1 * p + y1 * q)]
end

function rand_pol_with_np(P)
    _, x = polynomial_ring(GF(65521), "x" => 1:ambient_dim(P))
    ll = [Int.(v) for v in vertices(P)]
    return sum([rand(1:65520) * prod(x .^ l) for l in ll])
end

P = convex_hull([[15,0], [15,25], [0,25], [0,0]])
f = rand_pol_with_np(P)
g = rand_pol_with_np(P)
eqns = eci_from_dynsys(f, g)
eqns1 = eci_from_dynsys_for_mfp(f, g)
Q = get_eliminant_polytope(eqns)
Q1 = mixed_fiber_polytope(eqns1)
