struct AddNaiveInnerMostBlockLevelMap <: Gridap.Fields.Map
end

function Gridap.Arrays.return_cache(
  k::AddNaiveInnerMostBlockLevelMap,
  a::Gridap.Fields.ArrayBlock{<:ArrayBlock})

  # 1. Find a single element touched of a.array
  ai=testitem(a.array)

  # 2. Generate cache of k applied to ai
  ci=Gridap.Arrays.return_cache(k,ai)

  # 3. Generate array of caches
  cache_array=Array{typeof(ci),length(size(a.array))}(undef,size(a.array))
  for i in eachindex(a.array)
    if a.touched[i]
      cache_array[i]=Gridap.Arrays.return_cache(k,a.array[i])
    end
  end

  # 4. Generate ArrayBlock with the output entries
  bi       = Gridap.Arrays.evaluate!(ci,k,ai)
  barray   = Array{typeof(bi),length(size(a.array))}(undef,size(a.array))
  btouched = copy(a.touched)
  b=Gridap.Fields.ArrayBlock(barray,btouched)

  b,cache_array
end

function Gridap.Arrays.evaluate!(
  cache,
  k::AddNaiveInnerMostBlockLevelMap,
  a::Gridap.Fields.ArrayBlock{<:ArrayBlock}) where{A}

  b,cache_array=cache
  Gridap.Helpers.@check size(cache_array) == size(a.array)
  Gridap.Helpers.@check b.touched == a.touched
  for i in eachindex(a.array)
    if a.touched[i]
      b.array[i]=Gridap.Arrays.evaluate!(cache_array[i],k,a.array[i])
    end
  end
  b
end

function Gridap.Arrays.return_cache(
  k::AddNaiveInnerMostBlockLevelMap,
  a::Gridap.Fields.ArrayBlock{T}) where T

  n = length(size(a.array))
  Gridap.Helpers.@check n==1 || n==2

  # Generate an ArrayBlock with the same rank/shape as "a"
  # but 1x1 blocks instead of scalar entries
  if n==1
    ba_array=Array{VectorBlock{T}}(undef,size(a.array))
  else
    ba_array=Array{MatrixBlock{T}}(undef,size(a.array))
  end
  ba_touched=copy(a.touched)
  ba=Gridap.Fields.ArrayBlock(ba_array,ba_touched)

  # Generate 1x1 Naive ArrayBlocks
  for i in eachindex(ba_array)
    if ba_touched[i]
      if n==1
        sb_array=Vector{T}(undef,1)
        sb_touched=Vector{Bool}(undef,1)
      else
        sb_array=Matrix{T}(undef,1,1)
        sb_touched=Matrix{Bool}(undef,1,1)
      end
      sb_touched[1]=true
      sb=Gridap.Fields.ArrayBlock(sb_array,sb_touched)
      ba_array[i]=sb
    end
  end
  Gridap.Fields.ArrayBlock(ba_array,ba_touched)
end

function Gridap.Arrays.evaluate!(
  cache,
  k::AddNaiveInnerMostBlockLevelMap,
  a::Gridap.Fields.ArrayBlock{A}) where{A}
  Gridap.Helpers.@check size(cache.array) == size(a.array)
  Gridap.Helpers.@check cache.touched == a.touched
  for i in eachindex(a.array)
    if a.touched[i]
      sb=cache.array[i]
      sb.array[1]=a.array[i]
    end
  end
  cache
end
