module Tests
  using Test
  @time @testset "StaticCondensationMapTests" begin include("StaticCondensationMapTests.jl") end
  @time @testset "SumFacetMapTests" begin include("SumFacetMapTests.jl") end
  @time @testset "AddNaiveInnerMostBlockLevelMapTests" begin include("AddNaiveInnerMostBlockLevelMapTests.jl") end
  @time @testset "Scalar2ArrayBlockMapTests" begin include("Scalar2ArrayBlockMapTests.jl") end
  @time @testset "FEAutodiffTests" begin include("FEAutodiffTests.jl") end
  @time @testset "DarcyRTHTests" begin include("DarcyRTHTests.jl") end
  @time @testset "DarcyHDGTests" begin include("DarcyHDGTests.jl") end
  @time @testset "LinearElasticityHDGTests" begin include("LinearElasticityHDGTests.jl") end
  @time @testset "MultiFieldLagrangeMultipliersTests" begin include("MultiFieldLagrangeMultipliersTests.jl") end
  @time @testset "LocalFEOperatorTests" begin include("LocalFEOperatorTests.jl") end
  @time @testset "PoissonHHOTests" begin include("PoissonHHOTests.jl") end

end # module
