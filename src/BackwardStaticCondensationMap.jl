struct ReblockInteriorDofsMap <: Gridap.Fields.Map
end


# TO-DO: add ional keyword args to control the kind of linear system/factorization
# (i.e., symmetric and positive definite, symmetric indefinite, general unsymmetric, etc.)
struct BackwardStaticCondensationMap{IFT <: AbstractVector{<:Int},
                                     BFT <: Union{AbstractVector{<:Integer},
                                                  AbstractVector{<:AbstractVector{<:Integer}}}} <: Gridap.Fields.Map
  static_condensation::StaticCondensationMap{IFT,BFT}
  reblock::ReblockInteriorDofsMap
  function BackwardStaticCondensationMap(interior_fields::AbstractVector{<:Integer},
                              boundary_fields::AbstractVector{<:Integer})
    IFT=typeof(interior_fields)
    BFT=typeof(boundary_fields)
    new{IFT,BFT}(StaticCondensationMap(interior_fields,boundary_fields),ReblockInteriorDofsMap())
  end
  function BackwardStaticCondensationMap(interior_fields::AbstractVector{<:Integer},
                              boundary_fields::AbstractVector{<:AbstractVector{<:Integer}})
    # Implementation of BackwardStaticCondensation for multi-field Lagrange multipliers pending
    Gridap.Helpers.@notimplemented
  end
end

function Gridap.Arrays.return_cache(k::BackwardStaticCondensationMap{IFT,BFT},
                                    A::Gridap.Fields.MatrixBlock{<:Matrix{T}},
                                    b::Gridap.Fields.VectorBlock{<:Vector{T}},
                                    x::AbstractVector{T}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}
    cache1=Gridap.Arrays.return_cache(k.static_condensation,A,b)
    _,interior_brs,boundary_brs,_,_,_,_,_,=cache1
    cache2=Gridap.Arrays.return_cache(k.reblock,
                                      interior_brs,
                                      boundary_brs,
                                      x,
                                      x)
    (cache1,cache2)
end

function Gridap.Arrays.evaluate!(cache,
  k::BackwardStaticCondensationMap{IFT,BFT},
  A::Gridap.Fields.MatrixBlock{<:Matrix{T}},
  b::Gridap.Fields.VectorBlock{<:Vector{T}},
  x::AbstractVector{T}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}

   cache1,cache2=cache
   kdensify,interior_brs,boundary_brs,A11t,A21t,A12t,A22t,b1t,b2t=cache1
   A11_matblk,A11_cache=A11t
   _set_matblk!(A11_matblk,
                A,
                k.static_condensation.interior_fields,
                k.static_condensation.interior_fields)
   A12_matblk,A12_cache=A12t
   _set_matblk!(A12_matblk,
                A,
                k.static_condensation.interior_fields,
                k.static_condensation.boundary_fields)
   b1_vecblk,b1_cache=b1t
   b2_vecblk,b2_cache=b2t
   _set_vecblk!(b1_vecblk,b,k.static_condensation.interior_fields)
   _set_vecblk!(b2_vecblk,b,k.static_condensation.boundary_fields)
   A11=Gridap.Arrays.evaluate!(A11_cache,kdensify,interior_brs,interior_brs,A11_matblk)
   A12=Gridap.Arrays.evaluate!(A12_cache,kdensify,interior_brs,boundary_brs,A12_matblk)
   b1=Gridap.Arrays.evaluate!(b1_cache,kdensify,interior_brs,b1_vecblk)
   b2=Gridap.Arrays.evaluate!(b2_cache,kdensify,boundary_brs,b2_vecblk)

   b2 = x

   # b1 = b1-A12*b2
   LinearAlgebra.BLAS.gemv!('N',-1.0,A12,b2,1.0,b1)

   # TO-DO: ipiv is allocated on each call to getrf! :-(
   # A11 = L11 U11
   LUA11,ipiv,info=LinearAlgebra.LAPACK.getrf!(A11)
   Gridap.Helpers.@check info==0

   # b1 = inv(A11)*b1
   LinearAlgebra.LAPACK.getrs!('N', LUA11, ipiv, b1)

   Gridap.Arrays.evaluate!(cache2,k.reblock,interior_brs,boundary_brs,b1,b2)
end


function Gridap.Arrays.return_cache(k::ReblockInteriorDofsMap,
                                    interior_brs::Vector{<:Integer},
                                    boundary_brs::Vector{<:Integer},
                                    vinterior::Vector{T},
                                    vboundary::Vector{T}) where {T}
  Gridap.Helpers.@check length(boundary_brs) == 1
  nblks=length(interior_brs)+length(boundary_brs)
  array=Vector{Vector{T}}(undef,nblks)
  touched=Vector{Bool}(undef,nblks)
  touched .= true
  for iblk=1:length(interior_brs)
    array[iblk]=Vector{T}(undef,interior_brs[iblk])
  end
  for (biblk,iblk) in enumerate(length(interior_brs)+1:length(interior_brs)+length(boundary_brs))
    array[iblk]=Vector{T}(undef,boundary_brs[biblk])
  end
  Gridap.Fields.ArrayBlock(array,touched)
end

function Gridap.Arrays.evaluate!(cache,
                                 k::ReblockInteriorDofsMap,
                                 interior_brs::Vector{<:Integer},
                                 boundary_brs::Vector{<:Integer},
                                 vinterior::Vector{T},
                                 vboundary::Vector{T}) where {T}
  Gridap.Helpers.@check length(boundary_brs) == 1
  nblks=length(interior_brs)+length(boundary_brs)
  current=1
  for iblk=1:length(interior_brs)
    cache.array[iblk] .= vinterior[current:current+interior_brs[iblk]-1]
    current=current+interior_brs[iblk]
  end
  current=1
  for (biblk,iblk) in enumerate(length(interior_brs)+1:length(interior_brs)+length(boundary_brs))
    cache.array[iblk] .= vboundary[current:current+boundary_brs[biblk]-1]
    current=current+boundary_brs[biblk]
  end
  cache
end
