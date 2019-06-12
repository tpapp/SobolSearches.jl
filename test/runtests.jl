using SobolSearches, Test

@testset "simple consistency checks" begin
    y = [0.2, 0.7]
    f(x) = sum(abs2, x .- y)
    s = sobol_search(f, [-2, -1.0], [3, 1.0], 10, 10000)
    @test s[1].value == f(s[1].position)
    @test s[1].position ≈ y atol = 0.02
    @test IndexStyle(typeof(s)) ≡ IndexLinear()
    @test size(s) == (10, )
    @test length(s) == 10
end
