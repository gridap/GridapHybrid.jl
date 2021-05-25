module Tests
  using Test
  @testset "DensifyInnerMostBlockLevelMap" begin include("DensifyInnerMostBlockLevelMapTests.jl") end
  @testset "StaticCondensationMap" begin include("StaticCondensationMapTests.jl") end
  @testset "IntegrationTests" begin include("IntegrationTests.jl") end
end # module
