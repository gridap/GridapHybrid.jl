
##

## Arrays

# TO-THINK: is this reasonable from an efficiency point of view?
function _restrict_cell_array_block_to_block(x,block)
  #@assert false
  lazy_map(i->i[block],x)
end

function _restrict_cell_array_block_to_block(
       x::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock},block)
    @assert length(x.values)==1
    Fill(x.values[1][block],length(x))
end

function _restrict_cell_array_block_to_block(x::Fill{<:Gridap.Fields.VectorBlock},block)
  Fill(x.value[block],length(x))
end

function _restrict_cell_array_block_to_block(
    x::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},block)
  x.args[block]
end


function _restrict_cell_array_block_to_block(
  x::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BroadcastingFieldOpMap}},block)
  args = map(xi->_restrict_cell_array_block_to_block(xi,block), x.args)
  lazy_map(x.maps.value, args...)
end

function _restrict_cell_array_block_to_block(
  x::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.Broadcasting}},block)
  args = map(xi->_restrict_cell_array_block_to_block(xi,block), x.args)
  lazy_map(x.maps.value, args...)
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

function Gridap.Arrays.lazy_map(::typeof(∇),
  a::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap},
                             <:Gridap.Fields.VectorBlock{<:Gridap.Fields.AffineMap}})
  _block_arrays(map(a->lazy_map(∇,a), a.args)...)
end

# I cannot implement this imization. We end up summing MatrixBlocks of different types!!!
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
