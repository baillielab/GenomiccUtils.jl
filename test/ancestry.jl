module TestAncestry

using Test
using GenomiccUtils
using CSV
using DataFrames
using MLJBase
using StatisticalMeasures

TESTDIR = joinpath(pkgdir(GenomiccUtils), "test")

@testset "Test processing input/output of scope" begin
    # Test formatting of input file for SCOPE
    input_file = joinpath(TESTDIR, "assets", "scope_inputs", "kgp.frq.strat")
    output_file, io = mktemp()
    GenomiccUtils.format_stratified_freqs(input_file, output_file)
    close(io)
    scope_freqs = CSV.read(output_file, DataFrame; delim="\t")
    @test Set(scope_freqs.CLST) == Set([1, 2, 3, 4, 5])
    @test scope_freqs.MAF isa Vector{Float64}
    @test all(scope_freqs.MAC .== "N")
    @test all(scope_freqs.NCHROBS .== "N")
    # Test reading scope estimates
    scope_ancestry_file = joinpath(TESTDIR, "assets", "scope_inputs", "scope_resultsQhat.txt")
    n_indiv = 150
    fam = DataFrame(FID=1:150, IID=1:150)
    Q = GenomiccUtils.read_scope_ancestry_estimates(n_indiv, scope_ancestry_file)
    @test size(Q) == (n_indiv, 5)
    GenomiccUtils.assign_scope_ancestry_estimates!(fam, Q; threshold=0.6)
    # Check first few lines
    @test fam.Superpopulation[1:6] == ["ADMIXED", "ADMIXED", "EUR", "ADMIXED", "AFR", "AFR"]
end

@testset "Test get models" begin
    xgboost = GenomiccUtils.get_XGBoost_model(nfolds=5)
    @test xgboost.resampling isa StratifiedCV
    @test xgboost.resampling.nfolds == 5
    @test xgboost.measure == MisclassificationRate()
    @test hasproperty(xgboost.model, :xg_boost_classifier)

    logistic_model = GenomiccUtils.get_logistic_classifier(nfolds=5)
    @test logistic_model.resampling isa StratifiedCV
    @test logistic_model.resampling.nfolds == 5
    @test logistic_model.measure == MisclassificationRate()
    @test hasproperty(logistic_model.model, :logistic_classifier)
end

@testset "Test ml_ancestry_estimation" begin
    tmpdir = mktempdir()
    output = joinpath(tmpdir, "ancestry.csv")
    genotypes_prefix = joinpath(TESTDIR, "assets", "genotypes", "kgp.merged.unrelated.shared")
    kgp_pedigree_file = joinpath(TESTDIR, "assets", "genotypes", "kgp.pedigrees.txt")
    estimate_ancestry(genotypes_prefix, kgp_pedigree_file;
        mode="ml",
        output=output, 
        threshold=0.5,
        model="xgboost",
        npcs=5,
        nfolds=3
    )
    ancestries = CSV.read(output, DataFrame)
    @test size(ancestries) == (33, 8)
    @test names(ancestries) == ["FID", "IID", "Superpopulation", "AFR", "AMR", "EAS", "EUR", "SAS"]
    admixed = subset(ancestries, :Superpopulation => x -> x .== "ADMIXED")
    if size(admixed, 1) > 0
        @test maximum(Matrix(admixed[!, ["AFR", "AMR", "EAS", "EUR", "SAS"]])) < 0.5
    end
    for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]
        pop_df = subset(ancestries, :Superpopulation => x -> x .== pop)
        if size(admixed, 1) > 0
            @test all(pop_df[!, pop] .> 0.5)
        end
    end
end

end

true