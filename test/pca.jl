module TestPCA

using Test
using GenomiccUtils

TESTDIR = joinpath(pkgdir(GenomiccUtils), "test")

@testset "Test plink2_pca_approx" begin
    tmpdir = mktempdir()
    output_prefix = joinpath(tmpdir, "pca_output")
    input_prefix = joinpath(TESTDIR, "assets", "genotypes", "kgp.merged.unrelated.shared")
    pcs = GenomiccUtils.plink2_pca_approx(input_prefix; npcs=5, output_prefix=output_prefix)
    @test size(pcs) == (103, 7)
    @test names(pcs) == ["#FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5"]
end

end

true