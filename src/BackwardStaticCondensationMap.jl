# TO-DO: add optional keyword args to control the kind of linear system/factorization
# (i.e., symmetric and positive definite, symmetric indefinite, general unsymmetric, etc.)
struct BackwardStaticCondensationMap{IFT <: AbstractVector{<:Int},
                                     BFT <: Union{AbstractVector{<:Integer},
                                                  AbstractVector{<:AbstractVector{<:Integer}}}} <: Gridap.Fields.Map
  static_condensation::StaticCondensationMap{IFT,BFT}
  function BackwardStaticCondensationMap(interior_fields::AbstractVector{<:Integer},
                              boundary_fields::AbstractVector{<:Integer})
    IFT=typeof(interior_fields)
    BFT=typeof(boundary_fields)
    new{IFT,BFT}(StaticCondensationMap(interior_fields,boundary_fields))
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
    Gridap.Arrays.return_cache(k.static_condensation,A,b)
end


function Gridap.Arrays.evaluate!(cache,
  k::BackwardStaticCondensationMap{IFT,BFT},
  A::Gridap.Fields.MatrixBlock{<:Matrix{T}},
  b::Gridap.Fields.VectorBlock{<:Vector{T}},
  x::AbstractVector{T}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}

   k,interior_bs,boundary_bs,A11t,A21t,A12t,A22t,b1t,b2t=cache
   A11_matblk,A11_cache=A11t
   #A21_matblk,A21_cache=A21t
   A12_matblk,A12_cache=A12t
   #A22_matblk,A22_cache=A22t
   b1_vecblk,b1_cache=b1t
   b2_vecblk,b2_cache=b2t
   A11=Gridap.Arrays.evaluate!(A11_cache,k,interior_bs,interior_bs,A11_matblk)
   #A21=Gridap.Arrays.evaluate!(A21_cache,k,boundary_bs,interior_bs,A21_matblk)
   A12=Gridap.Arrays.evaluate!(A12_cache,k,interior_bs,boundary_bs,A12_matblk)
   #A22=Gridap.Arrays.evaluate!(A22_cache,k,boundary_bs,interior_bs,A22_matblk)
   b1=Gridap.Arrays.evaluate!(b1_cache,k,interior_bs,b1_vecblk)
   b2=Gridap.Arrays.evaluate!(b2_cache,k,boundary_bs,b2_vecblk)

   b2 = x

   # TO-DO: ipiv is allocated on each call to getrf! :-(
   # A11 = L11 U11
   LUA11,ipiv,info=LinearAlgebra.LAPACK.getrf!(A11)
   Gridap.Helpers.@check info==0

   # b1 = inv(A11)*b1
   LinearAlgebra.LAPACK.getrs!('N', LUA11, ipiv, b1)

   # b1 = b1-A12*b2
   LinearAlgebra.BLAS.gemv!('N',-1.0,A12,b2,1.0,b1)

   (b1,b2)
end
