struct SumFacetsMap <: Gridap.Fields.Map end

function Gridap.Arrays.return_cache(
  k::SumFacetsMap,
  a::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check all(a.touched)
  m=Gridap.Fields.BroadcastingFieldOpMap(+)
  c=Gridap.Arrays.return_cache(m,a.array[1],a.array[2])
  v=Vector{typeof(c)}(undef,length(a.array)-1)
  v[1]=c
  res=Gridap.Arrays.evaluate!(v[1],m,a.array[1],a.array[2])
  for i=3:length(a.array)
    v[i-1]=Gridap.Arrays.return_cache(m,res,a.array[i])
    res=Gridap.Arrays.evaluate!(v[i-1],m,res,a.array[i])
  end
  v
end

function Gridap.Arrays.evaluate!(
  cache,
  k::SumFacetsMap,
  a::Gridap.Fields.VectorBlock{A}) where{A}
  m=Gridap.Fields.BroadcastingFieldOpMap(+)
  Gridap.Helpers.@check all(a.touched)
  res=Gridap.Arrays.evaluate!(cache[1],m,a.array[1],a.array[2])
  for i=3:length(a.array)
    res=Gridap.Arrays.evaluate!(cache[i-1],m,res,a.array[i])
  end
  res
end
