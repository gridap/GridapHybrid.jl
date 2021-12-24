# TO-DO: add ional keyword args to control the kind of linear system/factorization
# (i.e., symmetric and positive definite, symmetric indefinite, general unsymmetric, etc.)
struct StaticCondensationMap{IFT <: AbstractVector{<:Int},
                             BFT <: Union{AbstractVector{<:Integer},
                                          AbstractVector{<:AbstractVector{<:Integer}}}} <: Gridap.Fields.Map
  interior_fields :: IFT
  boundary_fields :: BFT
  function StaticCondensationMap(interior_fields::AbstractVector{<:Integer},
                              boundary_fields::AbstractVector{<:Integer})
    Gridap.Helpers.@check _check_preconditions(interior_fields,boundary_fields)
    IFT=typeof(interior_fields)
    BFT=typeof(boundary_fields)
    new{IFT,BFT}(interior_fields,boundary_fields)
  end
  function StaticCondensationMap(interior_fields::AbstractVector{<:Integer},
                              boundary_fields::AbstractVector{<:AbstractVector{<:Integer}})
    # Implementation of StaticCondensation for multi-field Lagrange multipliers pending
    Gridap.Helpers.@notimplemented
  end
end

function _check_preconditions(interior_fields::AbstractVector{<:Integer},
                              boundary_fields::AbstractVector{<:Integer})
  num_fields=length(interior_fields)+length(boundary_fields)
  field_within_range(x) = x in 1:num_fields
  result   = true
  result   = result && all(field_within_range, interior_fields)
  result   = result && all(field_within_range, boundary_fields)
  touched  = Vector{Bool}(undef,num_fields)
  touched .= false
  for f in interior_fields
    result = result && !touched[f]
    touched[f]=true
  end
  for f in boundary_fields
    result = result && !touched[f]
    touched[f]=true
  end
  result = result & all(touched)
end

function Gridap.Arrays.return_cache(k::StaticCondensationMap{IFT,BFT},
  Ab::Tuple{<:Gridap.Fields.MatrixBlock{<:Matrix{T}},<:Gridap.Fields.VectorBlock{<:Vector{T}}}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}
  Gridap.Arrays.return_cache(k,Ab[1],Ab[2])
end

function Gridap.Arrays.return_cache(k::StaticCondensationMap{IFT,BFT},
                                    A::Gridap.Fields.MatrixBlock{<:Matrix{T}},
                                    b::Gridap.Fields.VectorBlock{<:Vector{T}}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}

  A11_matblk = _build_matblk(A,k.interior_fields,k.interior_fields)
  A21_matblk = _build_matblk(A,k.boundary_fields,k.interior_fields)
  A12_matblk = _build_matblk(A,k.interior_fields,k.boundary_fields)
  A22_matblk = _build_matblk(A,k.boundary_fields,k.boundary_fields)
  b1_vecblk  = _build_vecblk(b,k.interior_fields)
  b2_vecblk  = _build_vecblk(b,k.boundary_fields)
  brs, bcs   = _compute_brs_bcs(A)
  Gridap.Helpers.@check brs == bcs
  interior_brs = brs[k.interior_fields]
  boundary_brs = brs[k.boundary_fields]
  k=DensifyInnerMostBlockLevelMap()
  A11_cache=Gridap.Arrays.return_cache(k,interior_brs,interior_brs,A11_matblk)
  A21_cache=Gridap.Arrays.return_cache(k,boundary_brs,interior_brs,A21_matblk)
  A12_cache=Gridap.Arrays.return_cache(k,interior_brs,boundary_brs,A12_matblk)
  A22_cache=Gridap.Arrays.return_cache(k,boundary_brs,boundary_brs,A22_matblk)
  b1_cache =Gridap.Arrays.return_cache(k,interior_brs,b1_vecblk)
  b2_cache =Gridap.Arrays.return_cache(k,boundary_brs,b2_vecblk)

  (k, interior_brs, boundary_brs,
   (A11_matblk,A11_cache),
   (A21_matblk,A21_cache),
   (A12_matblk,A12_cache),
   (A22_matblk,A22_cache),
   (b1_vecblk,b1_cache),
   (b2_vecblk,b2_cache))
end

function _compute_brs_bcs(a::Gridap.Fields.MatrixBlock{<:Matrix})
  brs=Vector{Int}(undef, size(a)[1])
  bcs=Vector{Int}(undef, size(a)[2])
  for j=1:size(a)[2]
    for i=1:size(a)[1]
      if (a.touched[i,j])
        brs[i]=size(a.array[i,j],1)
        bcs[j]=size(a.array[i,j],2)
      end
    end
  end
  brs,bcs
end


function _build_matblk(a::Gridap.Fields.MatrixBlock{<:Matrix{T}},
                       rf::AbstractVector{<:Integer},
                       cf::AbstractVector{<:Integer}) where {T}
  s = (length(rf),length(cf))
  array    = Array{Matrix{T}}(undef,s)
  touched  = Matrix{Bool}(undef,s)
  touched .= false
  for (J,BJ) in enumerate(cf)
    for (I,BI) in enumerate(rf)
      touched[I,J]=a.touched[BI,BJ]
      if (touched[I,J])
        array[I,J]=a.array[BI,BJ]
      end
    end
  end
  Gridap.Fields.ArrayBlock(array,touched)
end


function _set_matblk!(a::Gridap.Fields.MatrixBlock{<:Matrix{T}},
                      b::Gridap.Fields.MatrixBlock{<:Matrix{T}},
                      rf::AbstractVector{<:Integer},
                      cf::AbstractVector{<:Integer}) where {T}
  for (J,BJ) in enumerate(cf)
    for (I,BI) in enumerate(rf)
      Gridap.Helpers.@check a.touched[I,J] == b.touched[BI,BJ]
      if (b.touched[BI,BJ])
        a.array[I,J]=b.array[BI,BJ]
      end
    end
  end
end

function _build_vecblk(a::Gridap.Fields.VectorBlock{<:Vector{T}},
                       F::AbstractVector{<:Integer}) where {T}
  s = length(F)
  array    = Vector{Vector{T}}(undef,s)
  touched  = Vector{Bool}(undef,s)
  touched .= false
  for (I,BI) in enumerate(F)
    touched[I]=a.touched[BI]
    if (touched[I])
      array[I]=a.array[BI]
    end
  end
  Gridap.Fields.ArrayBlock(array,touched)
end

function _set_vecblk!(a::Gridap.Fields.VectorBlock{<:Vector{T}},
                      b::Gridap.Fields.VectorBlock{<:Vector{T}},
                      F::AbstractVector{<:Integer}) where {T}
  for (I,BI) in enumerate(F)
    Gridap.Helpers.@check a.touched[I] == b.touched[BI]
    if (a.touched[I])
      a.array[I]=b.array[BI]
    end
  end
end

function Gridap.Arrays.evaluate!(cache,
  k::StaticCondensationMap{IFT,BFT},
  Ab::Tuple{<:Gridap.Fields.MatrixBlock{<:Matrix{T}},<:Gridap.Fields.VectorBlock{<:Vector{T}}}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}
  Gridap.Arrays.evaluate!(cache,k,Ab[1],Ab[2])
end

function Gridap.Arrays.evaluate!(cache,
  k::StaticCondensationMap{IFT,BFT},
  A::Gridap.Fields.MatrixBlock{<:Matrix{T}},
  b::Gridap.Fields.VectorBlock{<:Vector{T}}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}

   kdensify,interior_brs,boundary_brs,A11t,A21t,A12t,A22t,b1t,b2t=cache
   A11_matblk,A11_cache=A11t
   _set_matblk!(A11_matblk,A,k.interior_fields,k.interior_fields)
   A21_matblk,A21_cache=A21t
   _set_matblk!(A21_matblk,A,k.boundary_fields,k.interior_fields)
   A12_matblk,A12_cache=A12t
   _set_matblk!(A12_matblk,A,k.interior_fields,k.boundary_fields)
   A22_matblk,A22_cache=A22t
   _set_matblk!(A22_matblk,A,k.boundary_fields,k.boundary_fields)
   b1_vecblk,b1_cache=b1t
   b2_vecblk,b2_cache=b2t
   _set_vecblk!(b1_vecblk,b,k.interior_fields)
   _set_vecblk!(b2_vecblk,b,k.boundary_fields)
   A11=Gridap.Arrays.evaluate!(A11_cache,kdensify,interior_brs,interior_brs,A11_matblk)
   A21=Gridap.Arrays.evaluate!(A21_cache,kdensify,boundary_brs,interior_brs,A21_matblk)
   A12=Gridap.Arrays.evaluate!(A12_cache,kdensify,interior_brs,boundary_brs,A12_matblk)
   A22=Gridap.Arrays.evaluate!(A22_cache,kdensify,boundary_brs,boundary_brs,A22_matblk)
   b1=Gridap.Arrays.evaluate!(b1_cache,kdensify,interior_brs,b1_vecblk)
   b2=Gridap.Arrays.evaluate!(b2_cache,kdensify,boundary_brs,b2_vecblk)

   # TO-DO: ipiv is allocated on each call to getrf! :-(
   # A11 = L11 U11
   LUA11,ipiv,info=LinearAlgebra.LAPACK.getrf!(A11)
   Gridap.Helpers.@check info==0

   # A12 = inv(A11)*A12
   LinearAlgebra.LAPACK.getrs!('N', LUA11, ipiv, A12)

   # A22 = A22 - A21*A12
   LinearAlgebra.BLAS.gemm!('N','N',-1.0,A21,A12,1.0,A22)

   # b1 = inv(A11)*b1
   LinearAlgebra.LAPACK.getrs!('N', LUA11, ipiv, b1)

   # b2 = b2-A21*b1
   LinearAlgebra.BLAS.gemv!('N',-1.0,A21,b1,1.0,b2)

   #
   (A22,b2)
end
