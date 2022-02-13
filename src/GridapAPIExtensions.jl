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

# I cannot implement this optimization. We end up summing MatrixBlocks of different types!!!
# function Gridap.Arrays.lazy_map(k::Gridap.Fields.Broadcasting,
#                                 a::Gridap.Arrays.LazyArray{<:Fill{<:DensifyInnerMostBlockLevelMap}},
#                                 b::Gridap.Arrays.LazyArray{<:Fill{<:DensifyInnerMostBlockLevelMap}})
#   a_arg=a.args[1]
#   b_arg=b.args[1]
#   lazy_map(a.maps.value,lazy_map(k,a_arg,b_arg))
# end
