module Tests
  using Test
  @testset "StaticCondensationMap" begin include("StaticCondensationMapTests.jl") end
  @testset "CellBoundaryTests" begin include("CellBoundaryTests.jl") end
  @testset "TempSkeletonTests" begin include("TempSkeletonTests.jl") end
  @testset "LowLevelDarcyRTHTests" begin include("LowLevelDarcyRTHTests.jl") end
  @testset "DarcyHDGTests" begin include("DarcyHDGTests.jl") end
end # module
