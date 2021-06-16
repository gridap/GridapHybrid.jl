module Tests
  using Test
  @testset "DensifyInnerMostBlockLevelMap" begin include("DensifyInnerMostBlockLevelMapTests.jl") end
  @testset "StaticCondensationMap" begin include("StaticCondensationMapTests.jl") end
  @testset "DarcyRTHTests" begin include("DarcyRTHTests.jl") end
end # module
