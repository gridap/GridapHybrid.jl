struct SkeletonCompressedVector{T,G<:Gridap.Geometry.FaceToCellGlue} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
  ctype_lface_pindex_to_value::Vector{Vector{Vector{T}}}
  glue::G
end

function _compressed_vector_from_glue(::Type{T}, glue) where T
  ctype_to_vector_block=
    Vector{Gridap.Fields.VectorBlock{T}}(undef,length(glue.ctype_to_lface_to_ftype))
  for ctype=1:length(glue.ctype_to_lface_to_ftype)
     num_facets=length(glue.ctype_to_lface_to_ftype[ctype])
     v=Vector{T}(undef,num_facets)
     t=Vector{Bool}(undef,num_facets)
     t.=true
     ctype_to_vector_block[ctype]=Gridap.Fields.ArrayBlock(v,t)
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block, glue.cell_to_ctype)
end

function Gridap.Arrays.array_cache(a::SkeletonCompressedVector{T}) where {T}
   _compressed_vector_from_glue(T,a.glue)
end

function Base.getindex(a::SkeletonCompressedVector,cell::Integer)
  c=array_cache(a)
  Gridap.Arrays.getindex!(c,a,cell)
end

function Gridap.Arrays.getindex!(cache,a::SkeletonCompressedVector,cell::Integer)
   vb=cache[cell]
   ctype=a.glue.cell_to_ctype[cell]
   for lface=1:length(vb)
    p = a.glue.cell_to_lface_to_pindex.ptrs[cell]-1
    pindex = a.glue.cell_to_lface_to_pindex.data[p+lface]
    vb.array[lface]=a.ctype_lface_pindex_to_value[ctype][lface][pindex]
   end
   vb
end

Base.size(a::SkeletonCompressedVector) = (length(a.glue.cell_to_ctype),)

Base.IndexStyle(::Type{<:SkeletonCompressedVector}) = IndexLinear()

function Gridap.Arrays.lazy_map(k::Gridap.Fields.LinearCombinationMap,
                  ::Type{T},b::SkeletonCompressedVector,c::Fill) where T
  d = Gridap.Arrays.CompressedArray([c.value,],Fill(1,length(c)))
  lazy_map(k,T,b,d)
end

function Gridap.Arrays.lazy_map(k::Gridap.Fields.LinearCombinationMap,
                                ::Type{T},
                                b::SkeletonCompressedVector,
                                c::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock}) where T
  if c.ptrs === b.glue.cell_to_ctype || c.ptrs == b.glue.cell_to_ctype
    ET=eltype(T)
    ctype_lface_pindex_to_r = Vector{Vector{Vector{ET}}}(undef,length(b.ctype_lface_pindex_to_value))
    for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
      lface_pindex_to_r = Vector{Vector{ET}}(undef,length(lface_pindex_to_value))
      for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
        pindex_to_r = Vector{ET}(undef,length(pindex_to_value))
        ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
        if ftype != Gridap.Arrays.UNSET
          for (pindex, value) in enumerate(pindex_to_value)
            pindex_to_r[pindex] = k(value,c.values[ctype][lface])
          end
        end
        lface_pindex_to_r[lface] = pindex_to_r
      end
      ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
    end
    SkeletonCompressedVector(ctype_lface_pindex_to_r,b.glue)
  else
    Gridap.Helpers.@notimplemented
  end
end

function Gridap.Arrays.lazy_map(
  k::typeof(Gridap.Arrays.evaluate),
  ::Type{T},
  a::Fill,
  b::SkeletonCompressedVector{<:AbstractArray{<:Point}}) where T
  ET=eltype(T)
  Gridap.Helpers.@check length(a) == length(b)
  ctype_lface_pindex_to_r = Vector{Vector{Vector{ET}}}(undef,length(b.ctype_lface_pindex_to_value))
  for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
    lface_pindex_to_r = Vector{Vector{ET}}(undef,length(lface_pindex_to_value))
    for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
      pindex_to_r = Vector{ET}(undef,length(pindex_to_value))
      ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
      if ftype != Gridap.Arrays.UNSET
        for (pindex, value) in enumerate(pindex_to_value)
          pindex_to_r[pindex] = Gridap.Arrays.evaluate(a.value,value)
        end
      end
      lface_pindex_to_r[lface] = pindex_to_r
    end
    ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
  end
  SkeletonCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

function Gridap.Arrays.lazy_map(
  k::typeof(evaluate),
  ::Type{T},
  a::Gridap.Arrays.CompressedArray,
  b::SkeletonCompressedVector) where T

  Gridap.Helpers.@check length(a) == length(b)
  ET=eltype(T)
  ctype_lface_pindex_to_r = Vector{Vector{Vector{ET}}}(undef,length(b.ctype_lface_pindex_to_value))
  for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
    lface_pindex_to_r = Vector{Vector{ET}}(undef,length(lface_pindex_to_value))
    for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
      pindex_to_r = Vector{ET}(undef,length(pindex_to_value))
      ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
      if ftype != Gridap.Arrays.UNSET
        for (pindex, value) in enumerate(pindex_to_value)
          pindex_to_r[pindex] = evaluate(a.values[ctype][lface],value)
        end
      end
      lface_pindex_to_r[lface] = pindex_to_r
    end
    ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
  end
  SkeletonCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

function Gridap.Arrays.lazy_map(
  k::typeof(evaluate),
  ::Type{T},
  b::SkeletonCompressedVector,
  a::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock}) where T
  Gridap.Helpers.@check length(a) == length(b)
  ET=eltype(T)
  ctype_lface_pindex_to_r = Vector{Vector{Vector{ET}}}(undef,length(b.ctype_lface_pindex_to_value))
  for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
    lface_pindex_to_r = Vector{Vector{ET}}(undef,length(lface_pindex_to_value))
    for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
      pindex_to_r = Vector{ET}(undef,length(pindex_to_value))
      ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
      if ftype != Gridap.Arrays.UNSET
        for (pindex, value) in enumerate(pindex_to_value)
          pindex_to_r[pindex] = evaluate(value,a.values[ctype][lface])
        end
      end
      lface_pindex_to_r[lface] = pindex_to_r
    end
    ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
  end
  SkeletonCompressedVector(ctype_lface_pindex_to_r,b.glue)
end



struct SkeletonVectorFromFacetVector{T,G,C,V} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
  glue::G
  cell_wise_facets_ids::C
  facet_vector::V
  function SkeletonVectorFromFacetVector(glue,cell_wise_facets_ids,facet_vector)
    G=typeof(glue)
    C=typeof(cell_wise_facets_ids)
    V=typeof(facet_vector)
    T=eltype(V)
    new{T,G,C,V}(glue,cell_wise_facets_ids,facet_vector)
  end
end

Base.size(a::SkeletonVectorFromFacetVector) = (length(a.glue.cell_to_ctype),)

Base.IndexStyle(::Type{<:SkeletonVectorFromFacetVector}) = IndexLinear()

function Gridap.Arrays.array_cache(a::SkeletonVectorFromFacetVector{T}) where T
   fvc=array_cache(a.facet_vector)
   cwfc=array_cache(a.cell_wise_facets_ids)
   vbc=_compressed_vector_from_glue(T,a.glue)
   fvc=_compressed_vector_from_glue(typeof(fvc),a.glue)
   for ctype=1:length(a.glue.ctype_to_lface_to_ftype)
    num_facets=length(a.glue.ctype_to_lface_to_ftype[ctype])
    for lface=1:num_facets
      fvc.values[ctype][lface]=array_cache(a.facet_vector)
    end
  end
   fvc,cwfc,vbc
end

function Gridap.Arrays.getindex!(cache,a::SkeletonVectorFromFacetVector,cell::Integer)
  fvc,cwfc,vbc=cache
  vb=vbc[cell]
  cwf=getindex!(cwfc,a.cell_wise_facets_ids,cell)
  for (lfacet,gfacet) in enumerate(cwf)
    fv=getindex!(fvc[cell][lfacet],a.facet_vector,gfacet)
    vb.array[lfacet]=fv
  end
  vb
end

function Base.getindex(a::SkeletonVectorFromFacetVector,cell::Integer)
  c=array_cache(a)
  Gridap.Arrays.getindex!(c,a,cell)
end

struct SkeletonVectorFromSplitDoFsCellVector{T,A,B,C,D,E} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
  glue::A
  cell_wise_facets_ids::B
  cell_vector_dofs::C
  row_dofs_split::D
  col_dofs_split::E
  function SkeletonVectorFromSplitDoFsCellVector(glue,
                                                 cell_wise_facets_ids,
                                                 cell_vector_dofs,
                                                 row_dofs_split,
                                                 col_dofs_split)

    Gridap.Helpers.@check length(row_dofs_split)  == length(cell_wise_facets_ids[1])
    Gridap.Helpers.@check isa(col_dofs_split,UnitRange) || isa(col_dofs_split,AbstractVector{<:UnitRange}) 
    
    if isa(col_dofs_split,AbstractVector{<:UnitRange})
      Gridap.Helpers.@check length(col_dofs_split) == length(cell_wise_facets_ids[1])
    end 
    
    Gridap.Helpers.@check length(glue.cell_to_ctype) == length(cell_vector_dofs)

    T=typeof(view(cell_vector_dofs[1],1:1,1:1))
    A=typeof(glue)
    B=typeof(cell_wise_facets_ids)
    C=typeof(cell_vector_dofs)
    D=typeof(row_dofs_split)
    E=typeof(col_dofs_split)
    new{T,A,B,C,D,E}(glue,
                     cell_wise_facets_ids,
                     cell_vector_dofs,
                     row_dofs_split,
                     col_dofs_split)
  end
end

Base.size(a::SkeletonVectorFromSplitDoFsCellVector) = (length(a.glue.cell_to_ctype),)

Base.IndexStyle(::Type{<:SkeletonVectorFromSplitDoFsCellVector}) = IndexLinear()

function _generate_skel_vec_from_split_dofs_cell_vec_cache(a::SkeletonVectorFromSplitDoFsCellVector{T}) where T

  ctype_to_vector_block=
       Vector{Gridap.Fields.VectorBlock{T}}(undef,length(a.glue.ctype_to_lface_to_ftype))

  for ctype=1:length(a.glue.ctype_to_lface_to_ftype)
     cell=findfirst(x->x==ctype,a.glue.cell_to_ctype)
     num_facets=length(a.glue.ctype_to_lface_to_ftype[ctype])
     vf=Vector{T}(undef,num_facets)
     tf=Vector{Bool}(undef,num_facets)
     tf.=true
     ctype_to_vector_block[ctype]=Gridap.Fields.ArrayBlock(vf,tf)
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block,a.glue.cell_to_ctype)
end

function Gridap.Arrays.array_cache(a::SkeletonVectorFromSplitDoFsCellVector{T}) where T
  cvdc=array_cache(a.cell_vector_dofs)
  cwfc=array_cache(a.cell_wise_facets_ids)
  cea=_generate_skel_vec_from_split_dofs_cell_vec_cache(a)
  cvdc,cwfc,cea
end

function Gridap.Arrays.getindex!(cache,a::SkeletonVectorFromSplitDoFsCellVector,cell::Integer)
  cvdc,cwfc,cea=cache
  ctype=a.glue.cell_to_ctype[cell]
  result=cea.values[ctype]
  cwf=getindex!(cwfc,a.cell_wise_facets_ids,cell)
  dofs=getindex!(cvdc,a.cell_vector_dofs,cell)
  for (lfacet,gfacet) in enumerate(cwf)
    if (isa(a.col_dofs_split,UnitRange))
      fdofs=view(dofs,a.row_dofs_split[lfacet], a.col_dofs_split)
    else
      fdofs=view(dofs,a.row_dofs_split[lfacet], a.col_dofs_split[lfacet])
    end
    result.array[lfacet]=fdofs
  end
  result
end

function Base.getindex(a::SkeletonVectorFromSplitDoFsCellVector,cell::Integer)
  c=array_cache(a)
  Gridap.Arrays.getindex!(c,a,cell)
end

struct SkeletonOwnerNref{T,G,SF} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
    cell_boundary_nref::SkeletonCompressedVector{T,G}
    sign_flip::SF
end

function Gridap.Arrays.array_cache(a::SkeletonOwnerNref{T}) where {T}
  array_cache(a.cell_boundary_nref),array_cache(a.sign_flip)
end

function Base.getindex(a::SkeletonOwnerNref,cell::Integer)
 c=array_cache(a)
 Gridap.Arrays.getindex!(c,a,cell)
end

function Gridap.Arrays.getindex!(cache,a::SkeletonOwnerNref,cell::Integer)
  cnref,csf=cache
  nref=getindex!(cnref,a.cell_boundary_nref,cell)
  sf=getindex!(csf,a.sign_flip,cell)
  for i in eachindex(nref.array)
    if sf[i]
      nref.array[i]=-nref.array[i]
    end
  end
  nref
end

Base.size(a::SkeletonOwnerNref) = size(a.cell_boundary_nref)
Base.IndexStyle(::Type{<:SkeletonOwnerNref}) = IndexLinear()

function Gridap.Arrays.lazy_map(
  k::typeof(evaluate),
  ::Type{T},
  b::SkeletonVectorFromFacetVector,
  a::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock}) where T
  Gridap.Helpers.@check length(a) == length(b)
  af=_cell_lfacet_vector_to_facet_vector(b.glue,a)
  bf=b.facet_vector
  bfx=lazy_map(evaluate,bf,af)
  SkeletonVectorFromFacetVector(b.glue,b.cell_wise_facets_ids,bfx)
end

function _cell_lfacet_vector_to_facet_vector(
      glue,
      cell_lface::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock{T}}) where T
  nftypes=_count_ftypes(glue,cell_lface)
  values=Vector{T}(undef,nftypes)
  for ctype in eachindex(cell_lface.values)
    for lface=1:length(cell_lface.values[ctype])
      ftype = glue.ctype_to_lface_to_ftype[ctype][lface]
      values[ftype]=cell_lface.values[ctype][lface]
    end
  end
  Gridap.Arrays.CompressedArray(values,glue.face_to_ftype)
end

function _count_ftypes(glue,cell_lface::Gridap.Arrays.CompressedArray)
  touched=Dict{Int,Bool}()
  c=0
  for ctype in eachindex(cell_lface.values)
      for lface=1:length(cell_lface.values[ctype])
        ftype = glue.ctype_to_lface_to_ftype[ctype][lface]
        if ! haskey(touched,ftype)
            touched[ftype]=true
            c=c+1
         end
      end
  end
  c
end

function Gridap.Arrays.lazy_map(
  ::typeof(∇),
  b::SkeletonVectorFromFacetVector)
  ∇bf=lazy_map(∇,b.cell_vector)
  SkeletonVectorFromFacetVector(b.glue,b.cell_wise_facets_ids,∇bf)
end


struct SkeletonVectorFromCellVector{T,G,V} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
  glue::G
  cell_vector::V
  function SkeletonVectorFromCellVector(glue,cell_vector)
    G=typeof(glue)
    V=typeof(cell_vector)
    T=eltype(V)
    new{T,G,V}(glue,cell_vector)
  end
end

Base.size(a::SkeletonVectorFromCellVector) = (length(a.glue.cell_to_ctype),)

Base.IndexStyle(::Type{<:SkeletonVectorFromCellVector}) = IndexLinear()

function Gridap.Arrays.array_cache(a::SkeletonVectorFromCellVector{T}) where T
   cvc=array_cache(a.cell_vector)
   vbc=_compressed_vector_from_glue(T,a.glue)
   cvc,vbc
end

function Gridap.Arrays.getindex!(cache,a::SkeletonVectorFromCellVector,cell::Integer)
  cvc,vbc=cache
  vb=vbc[cell]
  fv=getindex!(cvc,a.cell_vector,cell)
  for lfacet=1:length(vb.touched)
    vb.array[lfacet]=fv
  end
  vb
end

function Base.getindex(a::SkeletonVectorFromCellVector,cell::Integer)
  c=array_cache(a)
  Gridap.Arrays.getindex!(c,a,cell)
end
