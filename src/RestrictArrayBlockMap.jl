struct RestrictArrayBlockMap{N} <: Gridap.Fields.Map
  blocks  :: Array{CartesianIndex{N},N}
  touched :: Array{Bool,N}
end

function RestrictArrayBlockMap(blocks::Vector{<:Int})
  RestrictArrayBlockMap(map(CartesianIndex,blocks),[true for i=1:length(blocks)])
end

function Gridap.Arrays.return_cache(k::RestrictArrayBlockMap{N},
                                    v::Gridap.Fields.ArrayBlock{T,N}) where {T,N}
  array=Array{T,N}(undef,size(k.blocks))
  touched=Array{Bool,N}(undef,size(k.blocks))
  for i in eachindex(k.blocks)
    if k.touched[i] && v.touched[i]
      touched[i] = true
    end
  end
  Gridap.Fields.ArrayBlock(array,touched)
end

function Gridap.Arrays.evaluate!(cache,
                                 k::RestrictArrayBlockMap{N},
                                 v::Gridap.Fields.ArrayBlock{T,N}) where {T,N}
  a=cache
  for i in eachindex(k.blocks)
    if a.touched[i]
      a.array[i]=v.array[k.blocks[i]]
    end
  end
  a
end
