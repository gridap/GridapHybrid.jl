struct StaticCondensationMap{IFT <: AbstractVector{<:Int},
                             BFT <: Union{AbstractVector{<:Integer},
                                          AbstractVector{<:AbstractVector{<:Integer}}}} <: Gridap.Fields.Map
  interior_fields :: IFT
  boundary_fields :: BFT
  function StaticCondensationMap(interior_fields::AbstractVector{<:Integer},
                              boundary_fields::AbstractVector{<:Integer})
    @assert _check_preconditions(interior_fields,boundary_fields)
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
                                    A::Gridap.Fields.MatrixBlock{<:Matrix{T}},
                                    b::Gridap.Fields.VectorBlock{<:Vector{T}}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}

  A11_matblk = _build_matblk(A,k.interior_fields,k.interior_fields)
  A21_matblk = _build_matblk(A,k.boundary_fields,k.interior_fields)
  A12_matblk = _build_matblk(A,k.interior_fields,k.boundary_fields)
  A22_matblk = _build_matblk(A,k.boundary_fields,k.boundary_fields)
  brs, bcs   = _compute_brs_bcs(A)
  @assert brs == bcs
  interior_brs = brs[k.interior_fields]
  boundary_brs = brs[k.boundary_fields]
  k=DensifyInnerMostBlockLevelMap()
  A11_cache=Gridap.Arrays.return_cache(k,interior_brs,interior_brs,A11_matblk)
  A21_cache=Gridap.Arrays.return_cache(k,boundary_brs,interior_brs,A21_matblk)
  A12_cache=Gridap.Arrays.return_cache(k,interior_brs,boundary_brs,A12_matblk)
  A22_cache=Gridap.Arrays.return_cache(k,boundary_brs,boundary_brs,A22_matblk)
  (k, interior_brs, boundary_brs,
   (A11_matblk,A11_cache),
   (A21_matblk,A21_cache),
   (A12_matblk,A12_cache),
   (A22_matblk,A22_cache))
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

function _build_vecblk(a::Gridap.Fields.VectorBlock{<:Vector{T}},
                       F::AbstractVector{<:Integer}) where {T}
  s = length(F)
  array    = Vector{Vector{T}}(undef,s)
  touched  = Vector{Bool}(undef,s)
  touched .= false
  for (I,BI) in enumerate(F)
    touched[I]=a.touched[BI]
    if (touched[I])
      array[I]=a.touched[BI]
    end
  end
  Gridap.VectorBlock(array,touched)
end


function Gridap.Arrays.evaluate!(cache,
  k::StaticCondensationMap{IFT,BFT},
  A::Gridap.Fields.MatrixBlock{<:Matrix{T}},
  b::Gridap.Fields.VectorBlock{<:Vector{T}}) where {IFT<:AbstractVector{<:Integer}, BFT <: AbstractVector{<:Integer}, T}

   k,interior_bs,boundary_bs,A11t,A21t,A12t,A22t=cache
   A11_matblk,A11_cache=A11t
   A21_matblk,A21_cache=A21t
   A12_matblk,A12_cache=A12t
   A22_matblk,A22_cache=A22t
   A11=Gridap.Arrays.evaluate!(A11_cache,k,interior_bs,interior_bs,A11_matblk)
   A21=Gridap.Arrays.evaluate!(A21_cache,k,boundary_bs,interior_bs,A21_matblk)
   A12=Gridap.Arrays.evaluate!(A12_cache,k,interior_bs,boundary_bs,A12_matblk)
   A22=Gridap.Arrays.evaluate!(A22_cache,k,boundary_bs,interior_bs,A22_matblk)
   (A11,A21,A12,A22)
end
