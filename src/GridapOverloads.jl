
##

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

# I cannot implement this optimization. We end up summing MatrixBlocks of different types!!!
# function Gridap.Arrays.lazy_map(k::Gridap.Fields.Broadcasting,
#                                 a::Gridap.Arrays.LazyArray{<:Fill{<:DensifyInnerMostBlockLevelMap}},
#                                 b::Gridap.Arrays.LazyArray{<:Fill{<:DensifyInnerMostBlockLevelMap}})
#   a_arg=a.args[1]
#   b_arg=b.args[1]
#   lazy_map(a.maps.value,lazy_map(k,a_arg,b_arg))
# end



## Polynomials

# Optimizing evaluation at a single point

function Gridap.Arrays.return_cache(f::Gridap.Polynomials.QCurlGradMonomialBasis{D,T},x::Point) where {D,T}
  ndof = size(f.qgrad)[1]
  r = Gridap.Arrays.CachedArray(zeros(VectorValue{D,T},(ndof,)))
  xs = [x]
  cf = Gridap.Arrays.return_cache(f,xs)
  r, cf, xs
end

function Gridap.Arrays.evaluate!(cache,f::Gridap.Polynomials.QCurlGradMonomialBasis{D,T},x::Point) where {D,T}
  r, cf, xs = cache
  xs[1] = x
  v = Gridap.Arrays.evaluate!(cf,f,xs)
  ndof = size(v,2)
  Gridap.Arrays.setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
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
