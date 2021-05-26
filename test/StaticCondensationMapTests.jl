module StaticCondensationMapTests
   using Test
   using Gridap
   using ExploringGridapHybridization

   x=rand(3,3)
   y=Array{Matrix{Float64}}(undef,(3,3))
   y[1,1]=x
   y[2,1]=x.+1
   y[3,1]=x.+2
   y[1,2]=x.+3
   y[1,3]=x.+5
   touched  = Array{Bool,2}(undef,(3,3))
   touched .= true
   touched[2:3,2:3].=false

   y_v=Vector{Vector{Float64}}(undef,3)
   touched_v=Vector{Bool}(undef,3)
   touched_v.=true
   x=rand(3)
   y_v[1]=x
   y_v[2]=x.+1
   y_v[3]=x.+2

   k=StaticCondensationMap([1,2],[3])
   cache=Gridap.Arrays.return_cache(k,
                                   Gridap.Fields.ArrayBlock(y,touched),
                                   Gridap.Fields.ArrayBlock(y_v,touched_v))

   a=Gridap.Arrays.evaluate!(cache,
                           k,
                           Gridap.Fields.ArrayBlock(y,touched),
                           Gridap.Fields.ArrayBlock(y_v,touched_v))
end
