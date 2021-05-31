struct DensifyInnerMostBlockLevelMap <: Gridap.Fields.Map
end

function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevelMap,
                                    a::Gridap.Fields.VectorBlock{<:Vector{T}}) where {T}
  Gridap.Helpers.@check all(a.touched)
  # Precompute the size of each row block.
  # We are assuming that the size of each row
  # block is going to be equivalent for all cells
  brs=Vector{Int}(undef,size(a)[1])
  for i=1:size(a)[1]
    s=s+length(a.array[i])
    brs[i]=length(a.array[i])
  end
  Gridap.Arrays.CachedArray(Vector{T}(undef,(s,))),brs
end

function Gridap.Arrays.return_cache(k   :: DensifyInnerMostBlockLevelMap,
                                    brs :: Vector{<:Integer},
                                    a   :: Gridap.Fields.VectorBlock{<:Vector{T}}) where {T}
  Gridap.Helpers.@check length(brs)==size(a)[1]
  Gridap.Arrays.CachedArray(Vector{T}(undef,(sum(brs),)))
end


function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevelMap,
                                    a::Gridap.Fields.MatrixBlock{<:Vector{T}}) where {T}
  Gridap.Helpers.@check _check_preconditions(k,a)
  s=Vector{Int}(undef,2)
  s[1]=0
  s[2]=size(a)[2]
  # Precompute the size of each row block.
  # We are assuming that the size of each row
  # block is going to be equivalent for all cells
  brs=Vector{Int}(undef,size(a)[1])
  for i=1:size(a)[1]
    for j=1:size(a)[2]
      if (a.touched[i,j])
        s[1]=s[1]+length(a.array[i,j])
        brs[i]=length(a.array[i,j])
        break
      end
    end
  end
  Gridap.Arrays.CachedArray(Array{T,2}(undef,Tuple(s))),brs
end

function Gridap.Arrays.return_cache(k   :: DensifyInnerMostBlockLevelMap,
                                    brs :: Vector{<:Integer},
                                    a   :: Gridap.Fields.MatrixBlock{<:Vector{T}}) where {T}
  Gridap.Helpers.@check length(brs)==size(a)[1]
  Gridap.Arrays.CachedArray(Array{T,2}(undef,(sum(brs),size(a)[2])))
end

function _check_preconditions(k::DensifyInnerMostBlockLevelMap,
                              a::Gridap.Fields.MatrixBlock{<:Vector{T}}) where {T}
  # Double check that, for each row, there is at least one non-zero block
  found=false
  for i=1:size(a)[1]
    found=false
    for j=1:size(a)[2]
      if (a.touched[i,j])
        found = true
        break
      end
    end
    if !found
      break
    end
  end
  found
end

function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevelMap,
                                    a::Gridap.Fields.VectorBlock{<:Matrix{T}}) where {T}
  Gridap.Helpers.@check all(a.touched)
  s=Vector{Int}(undef,2)
  s[1]=0
  s[2]=size(a.array[1])[2]
  for i=1:size(a)[1]
    s[1]=s[1]+size(a.array[i])[1]
  end
  Gridap.Arrays.CachedArray(Array{T,2}(undef,Tuple(s)))
end

function Gridap.Arrays.return_cache(k  ::DensifyInnerMostBlockLevelMap,
                                    brs::Vector{<:Integer},
                                    cs ::Integer,
                                    a  ::Gridap.Fields.VectorBlock{<:Matrix{T}}) where {T}
  Gridap.Arrays.CachedArray(Array{T,2}(undef,(sum(brs),cs)))
end


function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevelMap,
                                    a::Gridap.Fields.MatrixBlock{<:Matrix{T}}) where {T}
  Gridap.Helpers.@check _check_preconditions(k,a)
  s=Vector{Int}(undef,2)
  s[1]=0
  s[2]=0
  # Precompute the size of each row/col block.
  # We are assuming that the size of each row/col
  # block is going to be equivalent for all cells
  brs=Vector{Int}(undef,size(a)[1])
  bcs=Vector{Int}(undef,size(a)[2])
  # Row traversal
  for i=1:size(a)[1]
    for j=1:size(a)[2]
      if (a.touched[i,j])
        s[1]=s[1]+size(a.array[i,j])[1]
        brs[i]=size(a.array[i,j])[1]
        break
      end
    end
  end
  #Column traversal
  for j=1:size(a)[2]
     for i=1:size(a)[1]
       if (a.touched[i,j])
        s[2]=s[2]+size(a.array[i,j])[2]
        bcs[j]=size(a.array[i,j])[2]
        break
       end
     end
  end
  Gridap.Arrays.CachedArray(Array{T,2}(undef,Tuple(s))),brs,bcs
end

function _check_preconditions(k::DensifyInnerMostBlockLevelMap,
                              a::Gridap.Fields.MatrixBlock{<:Matrix{T}}) where {T}
  # Double check that, for each row,
  # there is at least one non-zero block
  found=false
  for i=1:size(a)[1]
    found=false
    for j=1:size(a)[2]
      if (a.touched[i,j])
        found = true
        break
      end
    end
    if !found
      break
    end
  end
  # Double check that, for each col,
  # there is at least one non-zero block
  if (found)
    for j=1:size(a)[2]
      found=false
      for i=1:size(a)[1]
        if (a.touched[i,j])
          found = true
          break
        end
      end
      if !found
        break
      end
    end
  end
  found
end

function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevelMap,
                                    brs::Vector{<:Integer},
                                    bcs::Vector{<:Integer},
                                    a::Gridap.Fields.MatrixBlock{<:Matrix{T}}) where {T}
  Gridap.Arrays.CachedArray(Array{T,2}(undef,(sum(brs),sum(bcs))))
end


# Most general version
function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevelMap,
                                    a::Gridap.Fields.ArrayBlock{<:Array{T,M},N}) where {T,M,N}
      Gridap.Helpers.@check all(a.touched)
      block_size=size(a)
      max_M_N=max(M,N)
      densified_size=Vector{Int}(undef,max_M_N)
      current_index=Vector{Int}(undef,N)
      for D=1:N
         current_index[1:D-1].=1
         current_index[D+1:N].=1
         densified_size[D]=0
         for I=1:block_size[D]
           current_index[D]=I
           if (D<=M)
             densified_size[D]+=size(a.array[current_index...])[D]
           else
             densified_size[D]+=1
           end
         end
      end
      if (M>N)
        current_index.=1
        for I=N+1:M
          densified_size[I]=size(a.array[current_index...])[I]
        end
      end
      Gridap.Arrays.CachedArray(Array{T,max_M_N}(undef,Tuple(densified_size)))
end

function Gridap.Arrays.evaluate!(cache,
                                 k::DensifyInnerMostBlockLevelMap,
                                 a::Gridap.Fields.VectorBlock{<:Vector{T}}) where {T}
  cache_array, brs = cache
  _evaluate!(cache_array,brs,a)
end

function Gridap.Arrays.evaluate!(cache,
                                 k   :: DensifyInnerMostBlockLevelMap,
                                 brs :: Vector{<:Integer},
                                 a   :: Gridap.Fields.VectorBlock{<:Vector{T}}) where {T}
  _evaluate!(cache,brs,a)
end


function Gridap.Arrays.evaluate!(cache,
                                 k::DensifyInnerMostBlockLevelMap,
                                 a::Gridap.Fields.MatrixBlock{<:Vector{T}}) where {T}
  cache_array, brs = cache
  _evaluate!(cache_array,brs,a)
end

function Gridap.Arrays.evaluate!(cache,
                                 k   :: DensifyInnerMostBlockLevelMap,
                                 brs :: Vector{<:Integer},
                                 a   :: Gridap.Fields.MatrixBlock{<:Vector{T}}) where {T}
  _evaluate!(cache,brs,a)
end

function _evaluate!(cache_array,brs,a)
  output  = cache_array.array
  output .= zero(eltype(output))
  current_i=1
  for i=1:size(a)[1]
     range = current_i:current_i+brs[i]-1
     if (a.touched[i])
        output[range] = a.array[i]
     end
     current_i = current_i + length(range)
  end
  output
end

function Gridap.Arrays.evaluate!(cache,
                                 k::DensifyInnerMostBlockLevelMap,
                                 a::Gridap.Fields.MatrixBlock{<:Matrix{T}}) where {T}
  cache_array, brs, bcs = cache
  _evaluate!(cache_array, brs, bcs,a)
end

function Gridap.Arrays.evaluate!(cache,
  k::DensifyInnerMostBlockLevelMap,
  brs::Vector{<:Integer},
  bcs::Vector{<:Integer},
  a::Gridap.Fields.MatrixBlock{<:Matrix{T}}) where {T}
  _evaluate!(cache, brs, bcs, a)
end

function _evaluate!(cache_array, brs, bcs, a)
  output  = cache_array.array
  output .= zero(eltype(output))
  current_j=1
  for j=1:size(a)[2]
    current_i=1
    range_j = current_j:current_j+bcs[j]-1
    for i=1:size(a)[1]
      range_i = current_i:current_i+brs[i]-1
      if (a.touched[i,j])
         output[range_i,range_j] = a.array[i,j]
      end
      current_i = current_i + brs[i]
    end
    current_j = current_j + bcs[j]
  end
  output
end

function Gridap.Arrays.evaluate!(cache,
                                 k::DensifyInnerMostBlockLevelMap,
                                 a::Gridap.Fields.VectorBlock{<:Matrix{T}}) where {T}
  Gridap.Helpers.@check all(a.touched)
  output = cache.array
  current_i=1
  n=size(a.array[1])[2]
  for i=1:size(a)[1]
    range=current_i:current_i+size(a.array[i])[1]-1
    output[range,1:n] = a.array[i]
    current_i=current_i+length(range)
  end
  output
end

# Most general version
function Gridap.Arrays.evaluate!(cache,
                                 k::DensifyInnerMostBlockLevelMap,
                                 a::Gridap.Fields.ArrayBlock{<:Array{T,M},N}) where {T,M,N}
      Gridap.Helpers.@check 1==0 "This function still BUGGY"
      Gridap.Helpers.@check all(a.touched)
      max_M_N=max(M,N)
      block_size=size(a)
      scalar_size=size(cache)
      output=cache.array
      current_block_ranges    = Vector{UnitRange{Int}}(undef,max_M_N)
      upper_left_entry_index  = Vector{Int}(undef,max_M_N)
      upper_left_entry_index .= 1
      cinds=CartesianIndices(size(a))
      for cind in cinds
        current_block_size=size(a.array[cind])
        for (p,s) in enumerate(current_block_size)
            current_block_ranges[p]=upper_left_entry_index[p]:upper_left_entry_index[p]+s-1
            upper_left_entry_index[p]=upper_left_entry_index[p]+s
            if (upper_left_entry_index[p]>scalar_size[p])
              upper_left_entry_index[p]=1
            end
        end
        if (M<N)
          for p=M+1:N
            current_block_ranges[p]=upper_left_entry_index[p]:upper_left_entry_index[p]
            upper_left_entry_index[p]=upper_left_entry_index[p]+1
            if (upper_left_entry_index[p]>scalar_size[p])
              upper_left_entry_index[p]=1
            end
          end
        end
        output[current_block_ranges...]=a.array[cind]
      end
      output
end

function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevelMap,
      a::Gridap.Fields.ArrayBlock{<:Gridap.Fields.ArrayBlock{T,M} where {T,M},N}) where {N}
    cache_touched=a.touched
    i=findfirst(isone, cache_touched)
    cache_block=Gridap.Arrays.return_cache(k,a.array[i])
    cache_array=Array{typeof(cache_block),N}(undef,size(a))
    output_array=Array{Gridap.Arrays.return_type(k,a.array[i]),N}(undef,size(a))
    linds=LinearIndices(size(a))
    cinds=CartesianIndices(size(a))
    while i != nothing
      cache_array[i]=cache_block
      if (linds[i]+1 <= length(linds))
        i=findnext(isone, cache_touched, cinds[linds[i]+1])
        if (i!=nothing)
          cache_block=Gridap.Arrays.return_cache(k,a.array[i])
        end
      else
        i=nothing
      end
    end
    Gridap.Fields.ArrayBlock(cache_array,cache_touched),
    Gridap.Fields.ArrayBlock(output_array,cache_touched),
    linds,
    cinds
end

function Gridap.Arrays.evaluate!(cache,
    k::DensifyInnerMostBlockLevelMap,
    a::Gridap.Fields.ArrayBlock{<:Gridap.Fields.ArrayBlock{T,M} where {T,M},N}) where {N}
    cache_array, output_array, linds, cinds = cache
    Gridap.Helpers.@check cache_array.touched == a.touched
    i=findfirst(isone, cache_array.touched)
    while i != nothing
      output_array.array[i]=Gridap.Arrays.evaluate!(cache_array.array[i],k,a.array[i])
      if (linds[i]+1 <= length(linds))
        i=findnext(isone, cache_array.touched, cinds[linds[i]+1])
      else
        i=nothing
      end
    end
    output_array
end
