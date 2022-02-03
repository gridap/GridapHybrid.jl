struct Scalar2ArrayBlockMap <: Gridap.Fields.Map
end

function Gridap.Arrays.return_cache(k::Scalar2ArrayBlockMap,
                                    Ab::Tuple{<:Matrix{T},<:Vector{T}},
                                    bs::AbstractVector{<:Integer}) where {T}
  Gridap.Arrays.return_cache(k,Ab[1],Ab[2],bs)
end


function Gridap.Arrays.return_cache(k::Scalar2ArrayBlockMap,
                                    A::Matrix{T},
                                    b::Vector{T},
                                    bs::AbstractVector{<:Integer}) where {T}
  nblocks=size(bs)[1]

  Ab   = Array{Matrix{T},2}(undef, (nblocks,nblocks) )
  tAb  = Array{Bool,2}(undef, (nblocks,nblocks) )
  tAb .= true
  for j=1:length(bs)
    for i=1:length(bs)
      Ab[i,j]=Array{T,2}(undef,bs[i],bs[j])
    end
  end

  bb   = Array{Vector{T},1}(undef, (nblocks,) )
  tbb  = Array{Bool,1}(undef, (nblocks,) )
  tbb .= true
  for i=1:length(bs)
    bb[i]=Array{T,1}(undef,bs[i])
  end
  Gridap.Fields.ArrayBlock(Ab,tAb), Gridap.Fields.ArrayBlock(bb,tbb)
end

function Gridap.Arrays.evaluate!(cache,
  k::Scalar2ArrayBlockMap,
  Ab::Tuple{<:Matrix{T},<:Vector{T}},
  bs::AbstractVector{<:Integer}) where {T}
  Gridap.Arrays.evaluate!(cache,k,Ab[1],Ab[2],bs)
end

function Gridap.Arrays.evaluate!(cache,
                                 k::Scalar2ArrayBlockMap,
                                 A::Matrix{T},
                                 b::Vector{T},
                                 bs::AbstractVector{<:Integer}) where {T}
  Ab,bb=cache
  startj=1
  for j=1:length(bs)
    endj=startj+bs[j]-1
    starti=1
    for i=1:length(bs)
      endi=starti+bs[i]-1
      v=view(A,starti:endi,startj:endj)
      Gridap.Helpers.@check size(Ab.array[i,j])==size(v)
      Ab.array[i,j] .= v
      starti=endi+1
    end
    startj=endj+1
  end

  starti=1
  for i=1:length(bs)
    endi=starti+bs[i]-1
    bb.array[i] .= view(b,starti:endi)
    starti=endi+1
  end
  Ab,bb
end
