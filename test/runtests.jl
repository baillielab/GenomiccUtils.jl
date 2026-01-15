using GenomiccUtils
using Test

TESTDIR = joinpath(pkgdir(GenomiccUtils), "test")

@testset "GenomiccUtils.jl" begin
    @test include(joinpath(TESTDIR, "ancestry.jl"))
    @test include(joinpath(TESTDIR, "pca.jl"))
end
