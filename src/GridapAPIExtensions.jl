## Arrays
function Gridap.Arrays.return_cache(
  ::typeof(evaluate),
  a::Gridap.Fields.VectorBlock,
  x::Gridap.Fields.VectorBlock)

   Gridap.Helpers.@check length(a.array) == length(x.array)
   Gridap.Helpers.@check a.touched == x.touched
   Gridap.Helpers.@check all(a.touched)
   T=Gridap.Arrays.return_type(evaluate,a.array[1],x.array[1])
   Tc=typeof(Gridap.Arrays.return_cache(evaluate,a.array[1],x.array[1]))
   r=Vector{T}(undef,length(a.array))
   rc=Vector{Tc}(undef,length(a.array))
   rc[1]=Gridap.Arrays.return_cache(evaluate,a.array[1],x.array[1])
   for i=2:length(a.array)
      if a.array[i] === a.array[1]
        rc[i]=rc[1]
      else
        rc[i]=Gridap.Arrays.return_cache(evaluate,a.array[i],x.array[i])
      end
   end
   (Gridap.Fields.ArrayBlock(r,a.touched),Gridap.Fields.ArrayBlock(rc,a.touched))
end

function Gridap.Arrays.evaluate!(
  cache,
  ::typeof(evaluate),
  a::Gridap.Fields.VectorBlock,
  x::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check length(a.array) == length(x.array)
  Gridap.Helpers.@check a.touched == x.touched
  Gridap.Helpers.@check all(a.touched)
  (r,rc)=cache
  Gridap.Helpers.@check length(r.array) == length(a.array)
  for i=1:length(a.array)
    r.array[i]=Gridap.Arrays.evaluate!(rc.array[i],evaluate,a.array[i],x.array[i])
  end
  r
end

function Gridap.Arrays.lazy_map(k::typeof(evaluate),
                                ::Type{T},
                                b::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock},
                                c::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock}) where T
  Gridap.Helpers.@check length(b) == length(c)
  if c.ptrs === b.ptrs || c.ptrs == b.ptrs
     values_r = Vector{T}(undef,length(b.values))
     for i=1:length(b.values)
        values_r[i] =evaluate(k,b.values[i],c.values[i])
     end
  else
    Gridap.Helpers.@notimplemented
  end
  Gridap.Arrays.CompressedArray(values_r,b.ptrs)
end

function Gridap.Arrays.lazy_map(k::typeof(evaluate),
                                ::Type{T},
                                b::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock},
                                c::Fill{<:Gridap.Fields.VectorBlock}) where T
  Gridap.Helpers.@check length(b) == length(c)
  values_r = Vector{T}(undef,length(b.values))
  for i=1:length(b.values)
    values_r[i]=evaluate(k,b.values[i],c.value)
  end
  Gridap.Arrays.CompressedArray(values_r,b.ptrs)
end

function Gridap.Arrays.lazy_map(
  ::typeof(evaluate), a::Gridap.Arrays.LazyArray{<:Fill{typeof(Gridap.Fields.linear_combination)}},
                      x::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock})

  i_to_values = a.args[1]
  i_to_basis = a.args[2]
  i_to_basis_x = lazy_map(evaluate,i_to_basis,x)
  lazy_map(Gridap.Fields.LinearCombinationMap(:),i_to_values,i_to_basis_x)
end

# The following three function overloads are required in order to be able to
# have a FE function in place of the trial operand in an integral over the
# cell boundaries
function Gridap.Fields.linear_combination(u::Vector,
                                          f::Gridap.Fields.ArrayBlock)
  i::Int = findfirst(f.touched)
  fi = f.array[i]
  ufi = Gridap.Fields.linear_combination(u,fi)
  g = Vector{typeof(ufi)}(undef,length(f.touched))
  for i in eachindex(f.touched)
    if f.touched[i]
      g[i] = Gridap.Fields.linear_combination(u,f.array[i])
    end
  end
  ArrayBlock(g,f.touched)
end

function Gridap.Arrays.return_cache(k::Gridap.Fields.LinearCombinationMap,
                                    u::Vector,
                                    fx::Gridap.Fields.ArrayBlock)
  i::Int = findfirst(fx.touched)
  fxi = fx.array[i]
  li = return_cache(k,u,fxi)
  ufxi = evaluate!(li,k,u,fxi)
  l = Vector{typeof(li)}(undef,size(fx.array))
  g = Vector{typeof(ufxi)}(undef,size(fx.array))
  for i in eachindex(fx.array)
    if fx.touched[i]
      l[i] = return_cache(k,u,fx.array[i])
    end
  end
  ArrayBlock(g,fx.touched),l
end

function Gridap.Arrays.evaluate!(cache,
                                 k::Gridap.Fields.LinearCombinationMap,
                                 u::Vector,
                                 fx::Gridap.Fields.ArrayBlock)
  g,l = cache
  Gridap.Helpers.@check g.touched == fx.touched
  for i in eachindex(fx.array)
    if fx.touched[i]
      g.array[i] = evaluate!(l[i],k,u,fx.array[i])
    end
  end
  g
end

function Gridap.Fields.linear_combination(u::Vector{<:Vector},
                                          f::Gridap.Fields.ArrayBlock)
  i::Int = findfirst(f.touched)
  ui = u[i]
  fi = f.array[i]
  uifi = Gridap.Fields.linear_combination(ui,fi)
  g = Vector{typeof(uifi)}(undef,length(f.touched))
  for i in eachindex(f.touched)
    if f.touched[i]
      g[i] = Gridap.Fields.linear_combination(u[i],f.array[i])
    end
  end
  ArrayBlock(g,f.touched)
end

function Gridap.Arrays.return_cache(k::Gridap.Fields.LinearCombinationMap,
                                    u::Vector{<:Vector},
                                    fx::Gridap.Fields.ArrayBlock)
  i::Int = findfirst(fx.touched)
  fxi = fx.array[i]
  li = return_cache(k,u[i],fxi)
  ufxi = evaluate!(li,k,u[i],fxi)
  l = Vector{typeof(li)}(undef,size(fx.array))
  g = Vector{typeof(ufxi)}(undef,size(fx.array))
  for i in eachindex(fx.array)
    if fx.touched[i]
      l[i] = return_cache(k,u[i],fx.array[i])
    end
  end
  ArrayBlock(g,fx.touched),l
end

function Gridap.Arrays.evaluate!(cache,
                                 k::Gridap.Fields.LinearCombinationMap,
                                 u::Vector,
                                 fx::Gridap.Fields.ArrayBlock)
  g,l = cache
  Gridap.Helpers.@check g.touched == fx.touched
  for i in eachindex(fx.array)
    if fx.touched[i]
      g.array[i] = evaluate!(l[i],k,u[i],fx.array[i])
    end
  end
  g
end

function Gridap.Arrays.return_cache(
  #::typeof(evaluate),
  a::Gridap.Fields.VectorBlock,
  x::Gridap.Fields.VectorBlock)

   Gridap.Helpers.@check length(a.array) == length(x.array)
   Gridap.Helpers.@check a.touched == x.touched
   Gridap.Helpers.@check all(a.touched)
   T=Gridap.Arrays.return_type(evaluate,a.array[1],x.array[1])
   Tc=typeof(Gridap.Arrays.return_cache(evaluate,a.array[1],x.array[1]))
   r=Vector{T}(undef,length(a.array))
   rc=Vector{Tc}(undef,length(a.array))
   rc[1]=Gridap.Arrays.return_cache(evaluate,a.array[1],x.array[1])
   for i=2:length(a.array)
      if a.array[i] === a.array[1]
        rc[i]=rc[1]
      else
        rc[i]=Gridap.Arrays.return_cache(evaluate,a.array[i],x.array[i])
      end
   end
   (Gridap.Fields.ArrayBlock(r,a.touched),Gridap.Fields.ArrayBlock(rc,a.touched))
end

function Gridap.Arrays.evaluate!(
  cache,
  #::typeof(evaluate),
  a::Gridap.Fields.VectorBlock,
  x::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check length(a.array) == length(x.array)
  Gridap.Helpers.@check a.touched == x.touched
  Gridap.Helpers.@check all(a.touched)
  (r,rc)=cache
  Gridap.Helpers.@check length(r.array) == length(a.array)
  for i=1:length(a.array)
    r.array[i]=Gridap.Arrays.evaluate!(rc.array[i],evaluate,a.array[i],x.array[i])
  end
  r
end


function Gridap.Arrays.return_cache(
  a::Gridap.Fields.Field,
  x::Gridap.Fields.VectorBlock)

   Gridap.Helpers.@check all(x.touched)
   T=Gridap.Arrays.return_type(evaluate,a,x.array[1])
   Tc=typeof(Gridap.Arrays.return_cache(evaluate,a,x.array[1]))
   r=Vector{T}(undef,length(x.array))
   rc=Vector{Tc}(undef,length(x.array))
   rc[1]=Gridap.Arrays.return_cache(evaluate,a,x.array[1])
   for i=2:length(x.array)
      if x.array[i] === x.array[1]
        rc[i]=rc[1]
      else
        rc[i]=Gridap.Arrays.return_cache(evaluate,a,x.array[i])
      end
   end
   (Gridap.Fields.ArrayBlock(r,x.touched),Gridap.Fields.ArrayBlock(rc,x.touched))
end

function Gridap.Arrays.evaluate!(
  cache,
  a::Gridap.Fields.Field,
  x::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check all(x.touched)
  (r,rc)=cache
  for i=1:length(x)
    r.array[i]=Gridap.Arrays.evaluate!(rc.array[i],evaluate,a,x.array[i])
  end
  r
end

function Gridap.Arrays.return_value(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.VectorBlock,
  h::Gridap.Fields.VectorBlock)
  evaluate(k,f,h)
end

function Gridap.Arrays.return_cache(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.VectorBlock,
  h::Gridap.Fields.VectorBlock) where {A,N}
  Gridap.Helpers.@check length(f)==length(h)
  Gridap.Helpers.@check f.touched==h.touched
  fi = Gridap.Arrays.testitem(f)
  hi = Gridap.Arrays.testitem(h)
  li = Gridap.Arrays.return_cache(k,fi,hi)
  fix = Gridap.Arrays.evaluate!(li,k,fi,hi)
  l = Vector{typeof(li)}(undef,size(h.array))
  g = Vector{typeof(fix)}(undef,size(h.array))
  for i in eachindex(h.array)
    if h.touched[i]
      l[i] = Gridap.Arrays.return_cache(k,f.array[i],h.array[i])
    end
  end
  Gridap.Fields.ArrayBlock(g,h.touched),l
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.VectorBlock,
  h::Gridap.Fields.VectorBlock)

  g,l = cache
  Gridap.Helpers.@check g.touched == h.touched
  for i in eachindex(h.array)
    if h.touched[i]
      g.array[i] = Gridap.Arrays.evaluate!(l[i],k,f.array[i],h.array[i])
    end
  end
  g
end

function Gridap.Arrays.return_value(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.Field,
  h::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}) where {N}
  evaluate(k,f,h)
end

function Gridap.Arrays.return_cache(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.Field,
  h::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}) where {N}
  hi = Gridap.Arrays.testitem(h)
  li = Gridap.Arrays.return_cache(k,f,hi)
  fix = Gridap.Arrays.evaluate!(li,k,f,hi)
  l = Array{typeof(li),N}(undef,size(h.array))
  g = Array{typeof(fix),N}(undef,size(h.array))
  for i in eachindex(h.array)
    if h.touched[i]
      l[i] = Gridap.Arrays.return_cache(k,f,h.array[i])
    end
  end
  Gridap.Fields.ArrayBlock(g,h.touched),l
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.Field,
  h::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}) where {N}
  g,l = cache
  Gridap.Helpers.@check g.touched == h.touched
  for i in eachindex(h.array)
    if h.touched[i]
      g.array[i] = Gridap.Arrays.evaluate!(l[i],k,f,h.array[i])
    end
  end
  g
end

function Gridap.Fields.constant_field(a::Gridap.Fields.ArrayBlock{T,N}) where {T<:Number,N}
  v=Array{Gridap.Fields.ConstantField{T},N}(undef,size(a))
  for (i,e) in enumerate(a.array)
    v[i]=Gridap.Fields.ConstantField(e)
  end
  Gridap.Fields.ArrayBlock(v,a.touched)
end

function Gridap.Arrays.return_value(
  op::Broadcasting{<:Operation},
  x::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}...) where N
  xi=map(a->a.array[1],x)
  T=Gridap.Arrays.return_type(op,xi...)
  v=Array{T,N}(undef,size(x[1]))
  for i in eachindex(v)
    xi=map(a->a.array[i],x)
    v[i]=Gridap.Arrays.return_value(op,xi...)
  end
  Gridap.Fields.ArrayBlock(v,x[1].touched)
end

function Gridap.Arrays.evaluate!(cache,
  op::Broadcasting{<:Operation},
  x::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}...) where N
  Gridap.Arrays.return_value(op,x...)
end

function Gridap.Arrays.return_cache(
  op::Broadcasting{<:Operation},
  x::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}...) where N
  xi=map(a->a.array[1],x)
  T=Gridap.Arrays.return_type(op,xi...)
  Tc=typeof(Gridap.Arrays.return_cache(op,xi...))
  v=Array{T,N}(undef,size(x[1]))
  vc=Array{Tc,N}(undef,size(x[1]))
  for i in eachindex(v)
    xi=map(a->a.array[i],x)
    vc[i]=Gridap.Arrays.return_cache(op,xi...)
  end
  Gridap.Fields.ArrayBlock(v,x[1].touched), Gridap.Fields.ArrayBlock(vc,x[1].touched)
end

function Base.:(∘)(a::Gridap.Fields.VectorBlock{<:Gridap.Fields.Field},
                b::Gridap.Fields.VectorBlock{<:Gridap.Fields.Field})
  Gridap.Helpers.@check size(a)==size(b)
  Gridap.Helpers.@check all(a.touched==b.touched)
  Gridap.Helpers.@check all(a.touched)
  T=Gridap.Arrays.return_type(∘,a.array[1],b.array[1])
  v=Vector{T}(undef,length(a.array))
  for i=1:length(v)
    v[i]=a.array[i]∘b.array[i]
  end
  Gridap.Fields.ArrayBlock(v,a.touched)
end

function Gridap.Arrays.return_cache(
  k::Gridap.Fields.IntegrationMap,
  aq::Gridap.Fields.VectorBlock{A},
  w::Gridap.Fields.VectorBlock{B},
  jq::Gridap.Fields.VectorBlock{C}) where{A,B,C}
  Gridap.Helpers.@check length(aq)==length(w)
  Gridap.Helpers.@check length(w)==length(jq)
  Gridap.Helpers.@check all(aq.touched)
  Gridap.Helpers.@check all(w.touched)
  Gridap.Helpers.@check all(jq.touched)
  aqi = Gridap.Arrays.testvalue(A)
  wi  = Gridap.Arrays.testvalue(B)
  jqi = Gridap.Arrays.testvalue(C)
  ci = Gridap.Arrays.return_cache(k,aqi,wi,jqi)
  hi = Gridap.Arrays.evaluate!(ci,k,aqi,wi,jqi)
  a = Vector{typeof(hi)}(undef,size(aq.array))
  b = Vector{typeof(ci)}(undef,size(aq.array))
  for i in eachindex(aq.array)
    b[i] = Gridap.Arrays.return_cache(k,
                        aq.array[i],
                        w.array[i],
                        jq.array[i])
  end
  Gridap.Fields.ArrayBlock(a,aq.touched), b
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Gridap.Fields.IntegrationMap,
  aq::Gridap.Fields.VectorBlock,
  w::Gridap.Fields.VectorBlock,
  jq::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check length(aq)==length(w)
  Gridap.Helpers.@check length(w)==length(jq)
  Gridap.Helpers.@check all(aq.touched)
  Gridap.Helpers.@check all(w.touched)
  Gridap.Helpers.@check all(jq.touched)
  a,b = cache
  # TO-DO: variable number of faces per cell boundary
  Gridap.Helpers.@check length(a)==length(aq)
  for i in eachindex(aq.array)
    a.array[i] = evaluate!(b[i],
                           k,
                           aq.array[i],
                           w.array[i],
                           jq.array[i])
  end
  a
end

function Gridap.Arrays.lazy_map(
  k   :: Gridap.Fields.IntegrationMap,
      :: Type{T},
  fx  :: AbstractArray{<:Gridap.Fields.ArrayBlock},
  w   :: AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractVector}},
  jtx :: AbstractArray{<:Gridap.Fields.ArrayBlock}) where T
  int=LazyArray(T,Fill(k,length(fx)),fx,w,jtx)
  sum_facets=lazy_map(SumFacetsMap(),int)
  lazy_map(DensifyInnerMostBlockLevelMap(),sum_facets)
end

function compute_bulk_to_skeleton_l2_projection_dofs(
  A::Matrix{<:Real},
  B::Matrix{<:Real})
  cache=return_cache(compute_bulk_to_skeleton_l2_projection_dofs,A,B)
  evaluate!(cache,compute_bulk_to_skeleton_l2_projection_dofs,A,B)
end

function Gridap.Arrays.return_cache(
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Matrix{<:Real},
  B::Matrix{<:Real})
  # c=CachedArray(B)
  # c
  nothing
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Matrix{<:Real},
  B::Matrix{<:Real})
  #setsize!(cache,size(B))
  #cache.array .= B
  #ldiv!(lu(A),cache.array) # TO-DO: explore in-place lu!()
  #cache.array
  A\B
end

function Gridap.Arrays.return_cache(
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Matrix{<:Real},
  B::Vector{<:Real})
  # c=CachedArray(B)
  # c
  nothing
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Matrix{<:Real},
  B::Vector{<:Real})
  #setsize!(cache,size(B))
  #cache.array .= B
  #ldiv!(lu(A),cache.array) # TO-DO: explore in-place lu!()
  #cache.array
  A\B
end

function setup_bulk_to_skeleton_l2_projected_fields(
  A::Array,
  B::AbstractArray{<:Gridap.Fields.Field})
  cache=return_cache(setup_bulk_to_skeleton_l2_projected_fields,A,B)
  evaluate!(cache,setup_bulk_to_skeleton_l2_projected_fields,A,B)
end

function Gridap.Arrays.return_cache(
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Array,
  B::AbstractArray{<:Gridap.Fields.Field})
  return_cache(Gridap.Fields.linear_combination,A,B)
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Array,
  B::AbstractArray{<:Gridap.Fields.Field})
  evaluate!(cache,Gridap.Fields.linear_combination,A,B)
end

function Gridap.Arrays.return_cache(
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Array,
  B::Transpose{T,<:AbstractArray{<:Gridap.Fields.Field}}) where T
  c=return_cache(k,A,B.parent)
  v=evaluate!(c,k,A,B.parent)
  return_cache(transpose,v), c
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Array,
  B::Transpose{T,<:AbstractArray{<:Gridap.Fields.Field}}) where T
  ct,clc=cache
  v=evaluate!(clc,Gridap.Fields.linear_combination,A,B.parent)
  evaluate!(ct,transpose,v)
end

# The following overloads are required to compute the L2
# projection of the bulk trial basis functions onto the space of traces
# (defined on the skeleton)

Tprojfuncs=Union{typeof(compute_bulk_to_skeleton_l2_projection_dofs),
        typeof(setup_bulk_to_skeleton_l2_projected_fields)}
function Gridap.Arrays.return_cache(
  k::Tprojfuncs,
  A::Gridap.Fields.VectorBlock,
  B::Gridap.Fields.VectorBlock)

  Gridap.Helpers.@check size(A.array) == size(B.array)
  Gridap.Helpers.@check A.touched == B.touched

  ai=testitem(A)
  xi=testitem(B)

  T=Gridap.Arrays.return_type(k,ai,xi)
  Tc=typeof(Gridap.Arrays.return_cache(k,ai,xi))
  r=Vector{T}(undef,size(A.array))
  rc=Vector{Tc}(undef,size(A.array))
  for i in eachindex(A)
    if A.touched[i]
      rc[i]=Gridap.Arrays.return_cache(k,A.array[i],B.array[i])
    end
  end
  (Gridap.Fields.ArrayBlock(r,A.touched),Gridap.Fields.ArrayBlock(rc,A.touched))
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Tprojfuncs,
  A::Gridap.Fields.VectorBlock,
  B::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check size(A.array) == size(B.array)
  Gridap.Helpers.@check A.touched == B.touched
  (r,rc)=cache
  Gridap.Helpers.@check size(r.array) == size(A.array)
  for i in eachindex(A)
    if A.touched[i]
      r.array[i]=Gridap.Arrays.evaluate!(rc.array[i],k,A.array[i],B.array[i])
    end
  end
  r
end


function Gridap.Arrays.return_cache(
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Vector,
  B::Gridap.Fields.MatrixBlock)
  nB=findall(B.touched)
  Gridap.Helpers.@check length(nB)==1
  xi=testitem(B)
  Gridap.Arrays.return_cache(k,A,xi)
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Vector,
  B::Gridap.Fields.MatrixBlock)
  nB=findall(B.touched)
  Gridap.Helpers.@check length(nB)==1
  xi=testitem(B)
  Gridap.Arrays.evaluate!(cache,k,A,xi)
end


# A [1,b1], e.g., [1,2] nonzero block (matrix to be inverted)
# B [1,b2], e.g., [1,3] nonzero block (several RHS)
function Gridap.Arrays.return_cache(
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Gridap.Fields.MatrixBlock,
  B::Gridap.Fields.MatrixBlock)

  Gridap.Helpers.@check size(A.array) == size(B.array)

  nA=findall(A.touched)
  nB=findall(B.touched)

  Gridap.Helpers.@check length(nA)==length(nB)
  Gridap.Helpers.@check length(nA)==1
  Gridap.Helpers.@check nA[1][1]==nB[1][1]

  tb=nB[1][2]
  nb=size(A.array)[1]
  ai=testitem(A)
  bi=testitem(B)
  T=Gridap.Arrays.return_type(k,ai,bi)
  touched  = Matrix{Bool}(undef,(1,nb))
  touched .= false
  touched[1,tb] = true
  r = Matrix{T}(undef,(1,nb))
  c = Gridap.Arrays.return_cache(k,ai,bi)
  Gridap.Fields.ArrayBlock(r,touched), c, tb, nb
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Gridap.Fields.MatrixBlock,
  B::Gridap.Fields.MatrixBlock)
  r,c,tb,nb=cache
  Gridap.Helpers.@check size(A.array) == size(B.array)
  Gridap.Helpers.@check size(A.array) == (nb,nb)
  Gridap.Helpers.@check size(r)       == (1,nb)
  Gridap.Helpers.@check length(findall(A.touched)) == 1
  Gridap.Helpers.@check length(findall(B.touched)) == 1
  Gridap.Helpers.@check r.touched[1,tb]
  ai=testitem(A)
  bi=testitem(B)
  r.array[1,tb]=evaluate!(c,k,ai,bi)
  r
end

# A [b,b], e.g.,  [8,8] nonzero block (matrix to be inverted)
# B [b  ], e.g.,  [8]   nonzero block (single RHS)
function Gridap.Arrays.return_cache(
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Gridap.Fields.MatrixBlock,
  B::Gridap.Fields.VectorBlock)

  Gridap.Helpers.@check size(A.array)[2] == size(B.array)[1]

  nA=findall(A.touched)
  nB=findall(B.touched)

  Gridap.Helpers.@check length(nA)==length(nB)
  Gridap.Helpers.@check length(nA)==1

  ai=testitem(A)
  bi=testitem(B)
  c=Gridap.Arrays.return_cache(k,ai,bi)
  c
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Gridap.Fields.MatrixBlock,
  B::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check length(findall(A.touched)) == 1
  Gridap.Helpers.@check length(findall(B.touched)) == 1
  ai=testitem(A)
  bi=testitem(B)
  evaluate!(cache,k,ai,bi)
end

# A [f,f]      , e.g., [2,2] nonzero block (matrix to be inverted)
# B [f]        , e.g., [2]   nonzero block (several RHS)
# Return [1]
function Gridap.Arrays.return_cache(
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Gridap.Fields.MatrixBlock{<:Matrix},
  B::Gridap.Fields.VectorBlock{<:Matrix})

  Gridap.Helpers.@check size(A.array)[1] == size(A.array)[2]
  Gridap.Helpers.@check size(A.array)[1] == size(A.array)[2]

  nA=findall(A.touched)
  nB=findall(B.touched)

  Gridap.Helpers.@check length(nA)==length(nB)
  Gridap.Helpers.@check length(nA)==1
  Gridap.Helpers.@check nA[1][1]==nB[1][1]
  Gridap.Helpers.@check nA[1][2]==nB[1][1]

  ai=testitem(A)
  bi=testitem(B)
  T=Gridap.Arrays.return_type(k,ai,bi)
  touched  = Vector{Bool}(undef,1)
  touched .= true
  r = Vector{T}(undef,1)
  c = Gridap.Arrays.return_cache(k,ai,bi)
  Gridap.Fields.ArrayBlock(r,touched), c
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(compute_bulk_to_skeleton_l2_projection_dofs),
  A::Gridap.Fields.MatrixBlock{<:Matrix},
  B::Gridap.Fields.VectorBlock{<:Matrix})
  r,c=cache
  Gridap.Helpers.@check size(A.array)[1] == size(A.array)[2]
  Gridap.Helpers.@check size(A.array)[1] == size(A.array)[2]
  Gridap.Helpers.@check length(findall(A.touched)) == 1
  Gridap.Helpers.@check length(findall(B.touched)) == 1
  Gridap.Helpers.@check r.touched[1]
  ai=testitem(A)
  bi=testitem(B)
  r.array[1]=evaluate!(c,k,ai,bi)
  r
end

function Gridap.Arrays.return_cache(
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Gridap.Fields.MatrixBlock,
  B::Gridap.Fields.MatrixBlock)

  Gridap.Helpers.@check size(A.array) == size(B.array)

  nA=findall(A.touched)
  nB=findall(B.touched)

  Gridap.Helpers.@check length(nA)==length(nB)
  Gridap.Helpers.@check length(nA)==1
  Gridap.Helpers.@check nA[1][1]==nB[1][1]

  tb=nA[1][2]
  nb=size(A.array)[2]
  ai=testitem(A)
  bi=testitem(B)
  T=Gridap.Arrays.return_type(k,ai,bi)
  touched  = Matrix{Bool}(undef,(1,nb))
  touched .= false
  touched[1,tb] = true
  r = Matrix{T}(undef,(1,nb))
  c = Gridap.Arrays.return_cache(k,ai,bi)
  Gridap.Fields.ArrayBlock(r,touched), c, tb, nb
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Gridap.Fields.MatrixBlock,
  B::Gridap.Fields.MatrixBlock)
  r,c,tb,nb=cache
  Gridap.Helpers.@check size(A.array) == size(B.array)
  Gridap.Helpers.@check size(r)       == (1,nb)
  Gridap.Helpers.@check length(findall(A.touched)) == 1
  Gridap.Helpers.@check length(findall(B.touched)) == 1
  Gridap.Helpers.@check r.touched[1,tb]
  ai=testitem(A)
  bi=testitem(B)
  r.array[1,tb]=evaluate!(c,k,ai,bi)
  r
end

# A: [1,f]
# B: [1]
# output [1]
function Gridap.Arrays.return_cache(
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Gridap.Fields.VectorBlock{<:Matrix},
  B::Gridap.Fields.MatrixBlock)

  Gridap.Helpers.@check size(A.array) == (1,)

  nB=findall(B.touched)
  Gridap.Helpers.@check length(nB)==1

  ai=testitem(A)
  bi=testitem(B)
  T=Gridap.Arrays.return_type(k,ai,bi)
  touched  = Vector{Bool}(undef,1)
  touched .= true
  r = Vector{T}(undef,1)
  c = Gridap.Arrays.return_cache(k,ai,bi)
  Gridap.Fields.ArrayBlock(r,touched), c
end

function Gridap.Arrays.evaluate!(
  cache,
  k::typeof(setup_bulk_to_skeleton_l2_projected_fields),
  A::Gridap.Fields.VectorBlock{<:Matrix},
  B::Gridap.Fields.MatrixBlock)
  r,c,=cache
  Gridap.Helpers.@check size(A.array) == (1,)
  Gridap.Helpers.@check length(findall(B.touched)) == 1
  Gridap.Helpers.@check length(findall(A.touched)) == 1
  Gridap.Helpers.@check r.touched[1]
  ai=testitem(A)
  bi=testitem(B)
  r.array[1]=evaluate!(c,k,ai,bi)
  r
end

# TO-DO: The current version of this function is very simple and
# it does not cover all possible scenarios. In particular, it only
# convers the case where all fields in uh are posed either on a cell
# triangulation or facet triangulation including resp. all cells and
# facets of the underlying background model.
function Gridap.FESpaces._compute_cell_ids(uh,ttrian::SkeletonTriangulation)
  collect(1:num_cells(ttrian))
end
