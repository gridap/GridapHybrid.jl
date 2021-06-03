
##

function Gridap.FESpaces.get_cell_shapefuns(model::DiscreteModel,
                       cell_reffe::AbstractArray{<:Gridap.ReferenceFEs.GenericRefFE{Gridap.ReferenceFEs.RaviartThomas}},
                       ::Gridap.ReferenceFEs.L2Conformity)
    Gridap.FESpaces.get_cell_shapefuns(model,cell_reffe,Gridap.ReferenceFEs.DivConformity())
end


## Arrays

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
