module AddNaiveInnerMostBlockLevelMapTests

  using Gridap
  using Gridap.Fields
  using ExploringGridapHybridization
  using Test

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
  afb       = ArrayBlock(vf,touchedf) # [f][b]

  afb1 = ExploringGridapHybridization.AddNaiveInnerMostBlockLevelMap()(afb) # [f][b][1]
  @test afb1[3][1][1] == afb[3][1]

  vf        = Vector{typeof(a)}(undef,4)
  touchedf  = Vector{Bool}(undef,4)
  touchedf .= true
  vf[1]     = a
  vf[2]     = a
  vf[3]     = a
  vf[4]     = a
  abf       = ArrayBlock(vf,touchedf) # [f]

  rfb1ama=Gridap.Fields.BroadcastingFieldOpMap(*)(afb1,abf) # [f][b][1]*[f] = [f][b][1]
  @test rfb1ama[1][1][1] == a .* a

  # First facet
  b=rand(2,4)
  v            = Matrix{typeof(a)}(undef,1,4)
  touched      = Matrix{Bool}(undef,1,4)
  touched     .= false
  touched[1,1] = true
  v[1,1]       = b
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
  v[1,2]       = b
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
  v[1,3]       = b
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
  v[1,4]       = b
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
  # [f][1,b][1,f]
  af1b1f    = ArrayBlock(vf,touchedf)

  # [f][b][1]*[f][1,b][1,f] = [f][b,b][1,f][1,f]
  res=Gridap.Fields.BroadcastingFieldOpMap(*)(rfb1ama,af1b1f)
  res[3][1,1][1,3] == (a .* a) .* b

  # First facet
  a=rand(2,4)
  v            = Vector{typeof(a)}(undef,4)
  touched      = Vector{Bool}(undef,4)
  touched     .= false
  touched[1]   = true
  v[1]         = a
  f1           = ArrayBlock(v,touched)

  v             = Vector{typeof(f1)}(undef,3)
  touched       = Vector{Bool}(undef,3)
  touched      .= false
  touched[1]    = true
  v[1,1]        = f1
  f1b           = ArrayBlock(v,touched)


  # Second facet
  v            = Vector{typeof(a)}(undef,4)
  touched      = Vector{Bool}(undef,4)
  touched     .= false
  touched[2]   = true
  v[2]         = a
  f2           = ArrayBlock(v,touched)

  v             = Vector{typeof(f2)}(undef,3)
  touched       = Vector{Bool}(undef,3)
  touched      .= false
  touched[1]    = true
  v[1]          = f2
  f2b           = ArrayBlock(v,touched)

  # Third facet
  v            = Vector{typeof(a)}(undef,4)
  touched      = Vector{Bool}(undef,4)
  touched     .= false
  touched[3]   = true
  v[3]         = a
  f3           = ArrayBlock(v,touched)

  v             = Vector{typeof(f3)}(undef,3)
  touched       = Vector{Bool}(undef,3)
  touched      .= false
  touched[1]    = true
  v[1]          = f3
  f3b           = ArrayBlock(v,touched)

  # Fourth facet
  v            = Vector{typeof(a)}(undef,4)
  touched      = Vector{Bool}(undef,4)
  touched     .= false
  touched[4]   = true
  v[4]         = a
  f4           = ArrayBlock(v,touched)

  v             = Vector{typeof(f4)}(undef,3)
  touched       = Vector{Bool}(undef,3)
  touched      .= false
  touched[1]    = true
  v[1]          = f4
  f4b           = ArrayBlock(v,touched)

  vf        = Vector{typeof(f1b)}(undef,4)
  touchedf  = Vector{Bool}(undef,4)
  touchedf .= true
  vf[1]     = f1b
  vf[2]     = f2b
  vf[3]     = f3b
  vf[4]     = f4b
  afbf      = ArrayBlock(vf,touchedf) # [f][b][f]


  b=rand(2,4)
  v            = Matrix{typeof(b)}(undef,1,3)
  touched      = Matrix{Bool}(undef,1,3)
  touched     .= false
  touched[1,1] = true
  v[1,1]       = b
  ab           = ArrayBlock(v,touched)

  vf        = Vector{typeof(ab)}(undef,4)
  touchedf  = Vector{Bool}(undef,4)
  touchedf .= true
  vf[1]     = ab
  vf[2]     = ab
  vf[3]     = ab
  vf[4]     = ab
  af1b       = ArrayBlock(vf,touchedf) # [f][1,b]

  af1b1 = ExploringGridapHybridization.AddNaiveInnerMostBlockLevelMap()(af1b) # [f][1,b][1]
  @test af1b1[3][1,1][1] == af1b[3][1,1]

  # [f][b][f]*[f][1,b][1] = [f][b,b][f]
  res=Gridap.Fields.BroadcastingFieldOpMap(*)(afbf,af1b1)
  @test res[2][1,1][2] == a .* b
end
