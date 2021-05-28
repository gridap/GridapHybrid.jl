
## Arrays

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

# TO-THINK: is this reasonable from an efficiency point of view?
function _restrict_cell_array_block_to_block(x,block)
  lazy_map(i->i[block],x)
end

# TO-THINK: Not sure if we are loosing generality by constraining a to be of type
#           LazyArray{...}. I need to restrict the cell-wise block array to each
#           individual block, and with a LazyArray{...} this is very efficient as
#           the array is already restricted to each block in the a.args member variable.
function Gridap.Arrays.lazy_map(::typeof(evaluate),
                                a::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
                                x::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractArray{<:Point}}})
  args = [ lazy_map(evaluate,a.args[pos],_restrict_cell_array_block_to_block(x,pos))
                                                           for pos=1:length(a.args) ]
  k = a.maps.value
  lazy_map(k,args...)
end

function Gridap.Arrays.lazy_map(::typeof(∇),
  a::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap},
                             <:Gridap.Fields.VectorBlock{<:Gridap.Fields.AffineMap}})
  _block_arrays(map(a->lazy_map(∇,a), a.args)...)
end

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

# function Gridap.Arrays.return_value(k::Gridap.Fields.Broadcasting{typeof(∘)},
#                                     g::Gridap.Fields.ArrayBlock{<:Transpose,N},
#                                     h::Gridap.Fields.Field) where N
#     i=findfirst(g.touched)
#     gi=g.array[i].parent
#     array=Array{typeof(gi),N}(undef,size(g.touched))
#     touched=g.touched
#     for i in findall(g.touched)
#       array[i]=g.array[i].parent
#     end
#     Transpose(Gridap.Arrays.return_value(k,Gridap.Fields.ArrayBlock(array,touched),h))
# end

@inline Gridap.Arrays.evaluate!(cache,f::Gridap.Fields.Broadcasting,
                 x::Transpose,y::Gridap.Fields.Field)=Transpose(evaluate(f,x.parent,y))
