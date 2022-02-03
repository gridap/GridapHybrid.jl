module Tests
  using Test
  @testset "StaticCondensationMapTests" begin include("StaticCondensationMapTests.jl") end
  @testset "SumFacetMapTests" begin include("SumFacetMapTests.jl") end
  @testset "AddNaiveInnerMostBlockLevelMapTests" begin include("AddNaiveInnerMostBlockLevelMapTests.jl") end
  @testset "Scalar2ArrayBlockMapTests" begin include("Scalar2ArrayBlockMapTests.jl") end
  @testset "CellBoundaryTests" begin include("CellBoundaryTests.jl") end
  @testset "DarcyRTHTests" begin include("DarcyRTHTests.jl") end
  @testset "DarcyHDGTests" begin include("DarcyHDGTests.jl") end
  @testset "MultiFieldLagrangeMultipliersTests" begin include("MultiFieldLagrangeMultipliersTests.jl") end
end # module
