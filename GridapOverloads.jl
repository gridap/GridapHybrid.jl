
## Arrays

# TO-DO: discuss with @fverdugo
# I needed to overload this function in order to be able to build the LazyArray below
#    lazy_map(Broadcasting(Reindex(np_array)),cell_wise_facets_ids)
function Gridap.Arrays.return_cache(f::Broadcasting,x::Union{Number,AbstractArray{<:Number}}...)
  s = map(Gridap.Arrays._size,x)
  bs = Base.Broadcast.broadcast_shape(s...)
  T = Gridap.Arrays.return_type(f.f,map(Gridap.Arrays.testitem,x)...)
  N = length(bs)
  r = fill(Gridap.Arrays.return_value(f.f,map(Gridap.Arrays.testitem,x)...),bs)
  cache = Gridap.Arrays.CachedArray(r)
  Gridap.Arrays._prepare_cache!(cache,x...)
  cache
end

# TO-THINK: is this reasonable from an efficiency point of view?
function _restrict_cell_array_block_to_block(x,block)
  lazy_map(i->i[block],x)
end

function Gridap.Arrays.lazy_map(::typeof(evaluate),
                                a::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
                                x::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractArray{<:Point}}})
  args = [ lazy_map(evaluate,a.args[pos],_restrict_cell_array_block_to_block(x,pos))
                                                           for pos=1:length(a.args) ]
  k = a.maps.value
  lazy_map(k,args...)
end

function Gridap.Arrays.lazy_map(::typeof(evaluate),::Type{T},a::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.PosNegReindex}},x::AbstractArray) where T
  i_to_iposneg = a.args[1]
  ipos_to_i = findall((x)->(x>0),i_to_iposneg)
  ineg_to_i = findall((x)->(x<0),i_to_iposneg)
  xpos = lazy_map(Reindex(x),ipos_to_i)
  xneg = lazy_map(Reindex(x),ineg_to_i)
  apos = lazy_map(Reindex(a.maps.value.values_pos),i_to_iposneg[ipos_to_i])
  aneg = lazy_map(Reindex(a.maps.value.values_neg),-i_to_iposneg[ineg_to_i])
  cpos = lazy_map(evaluate,apos,xpos)
  cneg = lazy_map(evaluate,aneg,xneg)
  function f(x)
    npos=0
    nneg=0
    o=similar(x)
    for (i,v)  in enumerate(x)
      if v>0
        npos+=1
        o[i]=npos
      else
        nneg+=1
        o[i]=-nneg
      end
    end
    o
  end
  lazy_map(Gridap.Fields.PosNegReindex(cpos,cneg),T,f(i_to_iposneg))
end
