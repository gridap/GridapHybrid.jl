module DensifyInnerMostBlockLevelMapTests
   using Test
   using Gridap
   using ExploringGridapHybridization

   x=rand(3)
   y=Array{Vector{Float64},2}(undef,(1,4))
   y[1,1]=x
   y[1,3]=x.+2
   y[1,4]=x.+3
   touched  = Array{Bool,2}(undef,(1,4))
   touched .= true
   touched[1,2]=false
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevelMap(),
                                    Gridap.Fields.ArrayBlock(y,touched))
   @test size(cache[1]) == (3,4)
   @test cache[2] == [3]
   Gridap.Arrays.evaluate!(cache,DensifyInnerMostBlockLevelMap(),
                           Gridap.Fields.ArrayBlock(y,touched))

   x=rand(1,3)
   y=Array{Matrix{Float64}}(undef,(4,))
   y[1]=x
   y[2]=x.+1
   y[3]=x.+2
   y[4]=x.+3
   touched  = Array{Bool,1}(undef,(4,))
   touched .= true
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevelMap(),
                                    Gridap.Fields.ArrayBlock(y,touched))
   Gridap.Arrays.evaluate!(cache,
                           DensifyInnerMostBlockLevelMap(),
                           Gridap.Fields.ArrayBlock(y,touched))

   @test size(cache) == (4,3)

   x=rand(2,3)
   y=Array{Matrix{Float64}}(undef,(3,3))
   y[1,1]=x
   y[2,1]=x.+1
   y[3,1]=x.+2
   y[1,2]=x.+3
   y[1,3]=x.+5
   touched  = Array{Bool,2}(undef,(3,3))
   touched .= true
   touched[2:3,2:3].=false
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevelMap(),
                                    Gridap.Fields.ArrayBlock(y,touched))

   y_v=Vector{Vector{Float64}}(undef,3)
   touched_v=Vector{Bool}(undef,3)
   x=rand(3)
   y_v[1]=x
   y_v[2]=x.+1
   y_v[3]=x.+2

   Gridap.Arrays.evaluate!(cache,
                           DensifyInnerMostBlockLevelMap(),
                           Gridap.Fields.ArrayBlock(y,touched))

   @test size(cache[1]) == (6,9)
   @test cache[2] == [2,2,2]
   @test cache[3] == [3,3,3]


   touched_parent      = Array{Bool,2}(undef,(3,3))
   touched_parent     .= false
   touched_parent[2,3] = true
   y_parent            = Array{Gridap.Fields.ArrayBlock{Matrix{Float64}},2}(undef,(3,3))
   y_parent[2,3]       = Gridap.Fields.ArrayBlock(y,touched)
   parent              = Gridap.Fields.ArrayBlock(y_parent,touched_parent)
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevelMap(),
                                    parent)
   result=Gridap.Arrays.evaluate!(cache,DensifyInnerMostBlockLevelMap(),
                              Gridap.Fields.ArrayBlock(y_parent,touched_parent))

end
