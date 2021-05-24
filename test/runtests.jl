module Tests
  using Test
  @testset "DensifyInnerMostBlockLevel" begin include("DensifyInnerMostBlockLevelTests.jl") end
  @testset "IntegrationTests" begin include("IntegrationTests.jl") end
end # module
