module SumFacetMapTests
   using Test
   using Gridap
   using GridapHybrid
   using Gridap.Fields
   using Gridap.Arrays

   m=GridapHybrid.SumFacetsMap()

   a=rand(2,4)
   v          = Vector{typeof(a)}(undef,3)
   touched    = Vector{Bool}(undef,3)
   touched   .= false
   touched[1] = true
   v[1]       = a
   ab         = ArrayBlock(v,touched)

   vf        = Vector{typeof(ab)}(undef,4)
   touchedf  = Vector{Bool}(undef,4)
   touchedf .= true
   vf[1]     = ab
   vf[2]     = ab
   vf[3]     = ab
   vf[4]     = ab
   abf       = ArrayBlock(vf,touchedf)

   c=return_cache(m,abf)
   result=evaluate!(c,m,abf)
   @test all(result[1]-4*a .≈ 0.0)


   # First facet
   a=rand(2,1)
   v            = Matrix{typeof(a)}(undef,1,4)
   touched      = Matrix{Bool}(undef,1,4)
   touched     .= false
   touched[1,1] = true
   v[1,1]       = a
   f1           = ArrayBlock(v,touched)

   v             = Matrix{typeof(f1)}(undef,1,3)
   touched       = Matrix{Bool}(undef,1,3)
   touched      .= false
   touched[1,1]  = true
   v[1,1]        = f1
   f1b           = ArrayBlock(v,touched)


   # Second facet
   v            = Matrix{typeof(a)}(undef,1,4)
   touched      = Matrix{Bool}(undef,1,4)
   touched     .= false
   touched[1,2] = true
   v[1,2]       = a
   f2           = ArrayBlock(v,touched)

   v             = Matrix{typeof(f2)}(undef,1,3)
   touched       = Matrix{Bool}(undef,1,3)
   touched      .= false
   touched[1,1]  = true
   v[1,1]        = f2
   f2b           = ArrayBlock(v,touched)


   # Third facet
   v            = Matrix{typeof(a)}(undef,1,4)
   touched      = Matrix{Bool}(undef,1,4)
   touched     .= false
   touched[1,3] = true
   v[1,3]       = a
   f3           = ArrayBlock(v,touched)

   v             = Matrix{typeof(f3)}(undef,1,3)
   touched       = Matrix{Bool}(undef,1,3)
   touched      .= false
   touched[1,1]  = true
   v[1,1]        = f3
   f3b           = ArrayBlock(v,touched)

   # Fourth facet
   v            = Matrix{typeof(a)}(undef,1,4)
   touched      = Matrix{Bool}(undef,1,4)
   touched     .= false
   touched[1,4] = true
   v[1,4]       = a
   f4           = ArrayBlock(v,touched)

   v             = Matrix{typeof(f4)}(undef,1,3)
   touched       = Matrix{Bool}(undef,1,3)
   touched      .= false
   touched[1,1]  = true
   v[1,1]        = f4
   f4b           = ArrayBlock(v,touched)

   vf        = Vector{typeof(f1b)}(undef,4)
   touchedf  = Vector{Bool}(undef,4)
   touchedf .= true
   vf[1]     = f1b
   vf[2]     = f2b
   vf[3]     = f3b
   vf[4]     = f4b
   abf       = ArrayBlock(vf,touchedf)

   m=GridapHybrid.SumFacetsMap()
   c=return_cache(m,abf)
   result=evaluate!(c,m,abf)
   c=return_cache(DensifyInnerMostBlockLevelMap(),result)
   result=evaluate!(c,DensifyInnerMostBlockLevelMap(),result)
   @test all(result[1,1]-hcat(a,a,a,a) .≈ 0.0)
end
