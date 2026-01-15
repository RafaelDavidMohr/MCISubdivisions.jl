using MCISubdivisions
using Oscar
using Test

const MCIS = MCISubdivisions

function random_linear_as_mci(n::Int)
    A = [i == j ? 1 : 0 for i in 1:n, j in 1:n]
    A = hcat(zeros(Int, n), A)
    V = (QQ).(rand(-100:100, n, n + 1))
    M = MCIS.MCI(V, A)
    return M
end

function hexagon_example()
    A = [0 1 2 2 1 0; 0 0 1 2 2 1]
    r = rand(1:100, 1, 2)
    g = rand(1:100, 1, 2)
    b = rand(1:100, 1, 2)
    V = hcat([r[1]],g,b,[r[2]])
    V = (QQ).(vcat(V, hcat([r[1]], 2*g, 3*b, [r[2]])))
    M = MCIS.MCI(V, A)
    return M
end

@testset "Setup" begin
    # random linear equations
    M = random_linear_as_mci(4)
    m = MCIS.MixedCell([collect(1:5)], M)
    wd = MCIS.WalkData(M, [m])
    @test isempty(wd.walls)

    # hexagon example
    M = hexagon_example()
    m = MCIS.MixedCell([[1,6], [2,3]], M)
    wd = MCIS.WalkData(M, [m])
    @test length(keys(wd.walls)) == 3
    @test length(findall(im -> first(im)[2] == 1, wd.walls)) == 1
    @test length(findall(im -> first(im)[2] == 2, wd.walls)) == 2
end

@testset "Mixed Cell Functions" begin
    # random linear equations
    M = random_linear_as_mci(4)
    d = (QQ).(rand(-100:100, 5))
    d0 = first(d)
    w = (QQ).([d0 - d[i] for i in 2:5])
    mixed_cell = MCIS.find_dual_tropical_root(M, d, w)
    @test mixed_cell.inds == [collect(1:5)]
    @test MCIS.is_partial_mixed_cell(M, mixed_cell)

    # hexagon example
    M = hexagon_example()
    inds = [1,6]
    for S in [[2,3], [4,5]]
        @test MCIS.is_partial_mixed_cell(M, MCIS.MixedCell([inds, S], M))
    end
    @test !MCIS.is_partial_mixed_cell(M, MCIS.MixedCell([inds, [3,4]], M))

    d = (QQ).([14,27,56,63,50,27])
    w = (QQ).([-16,-13])
    m = MCIS.find_dual_tropical_root(M, d, w)
    @test m.inds == [[1,6],[2,3]]
end

@testset "Mixed Cell Flips" begin
    # hexagon example
    M = hexagon_example()
    m = MCIS.MixedCell([[1,6], [2,3]], M)
    wd = MCIS.WalkData(M, [m])
    wls = collect(keys(wd.walls))
    c = wls[findfirst(c -> c.nzinds == [1,2,3,6], wls)]
    new_mc = MCIS.mixed_cell_flip(m, c, M, 1, false)
    @test length(new_mc) == 1
    @test first(new_mc).inds == [[2,3],[1,6]]
    m = MCIS.MixedCell([[1,6], [2,4]], M)
    wd = MCIS.WalkData(M, [m])
    wls = collect(keys(wd.walls))
    c = wls[findfirst(c -> c.nzinds == [1,3,4,6], wls)]
    new_mc = MCIS.mixed_cell_flip(m, c, M, 2, false)
    @test length(new_mc) == 1
    @test first(new_mc).inds == [[1,6],[2,3]]
    c = wls[findfirst(c -> c.nzinds == [1,2,4,6], wls)]
    new_mc = MCIS.mixed_cell_flip(m, c, M, 1, false)
    @test length(new_mc) == 1
    @test first(new_mc).inds == [[2,4,6]]
end

@testset "Mixed Volume" begin
    A = [0 1 2 2 1 0; 0 0 1 2 2 1]
    r = rand(1:100, 1, 2)
    g = rand(1:100, 1, 2)
    b = rand(1:100, 1, 2)
    V = hcat([r[1]],g,b,[r[2]])
    V = (QQ).(vcat(V, hcat([r[1]], 2*g, 3*b, [r[2]])))
    @test mixed_volume(A, V) == 3
end

#=
test case:
A = [5 4 3 2 1 0]
V = [1 1 1 1 1 1]
lift = [-60, -74, -26, 24, 31, -47]
=#
