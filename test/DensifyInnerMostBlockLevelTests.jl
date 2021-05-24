module DensifyInnerMostBlockLevelTests
   using Test
   using Gridap
   using ExploringGridapHybridization

   x=rand(3)
   y=Array{Vector{Float64},2}(undef,(1,4))
   y[1,1]=x
   y[1,2]=x.+1
   y[1,3]=x.+2
   y[1,4]=x.+3
   touched  = Array{Bool,2}(undef,(1,4))
   touched .= true
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevel(),
                                    Gridap.Fields.ArrayBlock(y,touched))
   @test size(cache) == (3,4)
   Gridap.Arrays.evaluate!(cache,DensifyInnerMostBlockLevel(),
                           Gridap.Fields.ArrayBlock(y,touched))

   x=rand(1,3)
   y=Array{Matrix{Float64}}(undef,(4,))
   y[1]=x
   y[2]=x.+1
   y[3]=x.+2
   y[4]=x.+3
   touched  = Array{Bool,1}(undef,(4,))
   touched .= true
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevel(),
                                    Gridap.Fields.ArrayBlock(y,touched))
   Gridap.Arrays.evaluate!(cache,
                           DensifyInnerMostBlockLevel(),
                           Gridap.Fields.ArrayBlock(y,touched))

   @test size(cache) == (4,3)

   x=rand(2,3)
   y=Array{Matrix{Float64}}(undef,(3,2))
   y[1,1]=x
   y[2,1]=x.+1
   y[3,1]=x.+2
   y[1,2]=x.+3
   y[2,2]=x.+4
   y[3,2]=x.+5
   touched  = Array{Bool,2}(undef,(3,2))
   touched .= true
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevel(),
                                    Gridap.Fields.ArrayBlock(y,touched))
   Gridap.Arrays.evaluate!(cache,
                           DensifyInnerMostBlockLevel(),
                           Gridap.Fields.ArrayBlock(y,touched))

  @test size(cache) == (6,6)

   touched_parent      = Array{Bool,2}(undef,(3,3))
   touched_parent     .= false
   touched_parent[2,3] = true
   y_parent            = Array{Gridap.Fields.ArrayBlock{Matrix{Float64}},2}(undef,(3,3))
   y_parent[2,3]       = Gridap.Fields.ArrayBlock(y,touched)
   parent              = Gridap.Fields.ArrayBlock(y_parent,touched_parent)
   cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevel(),
                                    parent)
   result=Gridap.Arrays.evaluate!(cache,DensifyInnerMostBlockLevel(),
                              Gridap.Fields.ArrayBlock(y_parent,touched_parent))

end
