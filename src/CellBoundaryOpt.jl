using Gridap
using FillArrays
using Gridap.Geometry


struct CellBoundaryCompressedVector{T,G<:Gridap.Geometry.FaceToCellGlue} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
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

function Gridap.Arrays.array_cache(a::CellBoundaryCompressedVector{T}) where {T}
   _compressed_vector_from_glue(T,a.glue)
end

function Base.getindex(a::CellBoundaryCompressedVector,cell::Integer)
  c=array_cache(a)
  Gridap.Arrays.getindex!(c,a,cell)
end

function Gridap.Arrays.getindex!(cache,a::CellBoundaryCompressedVector,cell::Integer)
   vb=cache[cell]
   ctype=a.glue.cell_to_ctype[cell]
   for lface=1:length(vb)
    p = a.glue.cell_to_lface_to_pindex.ptrs[cell]-1
    pindex = a.glue.cell_to_lface_to_pindex.data[p+lface]
    vb.array[lface]=a.ctype_lface_pindex_to_value[ctype][lface][pindex]
   end
   vb
end

Base.size(a::CellBoundaryCompressedVector) = (length(a.glue.cell_to_ctype),)

Base.IndexStyle(::Type{<:CellBoundaryCompressedVector}) = IndexLinear()

struct CellBoundaryVectorFromFacetVector{T,G,C,V} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
  eltype::Type{T}
  glue::G
  cell_wise_facets_ids::C
  facet_vector::V
  function CellBoundaryVectorFromFacetVector(glue,cell_wise_facets_ids,facet_vector)
    G=typeof(glue)
    C=typeof(cell_wise_facets_ids)
    V=typeof(facet_vector)
    T=eltype(V)
    new{T,G,C,V}(T,glue,cell_wise_facets_ids,facet_vector)
  end
end

Base.size(a::CellBoundaryVectorFromFacetVector) = (length(a.glue.cell_to_ctype),)

Base.IndexStyle(::Type{<:CellBoundaryVectorFromFacetVector}) = IndexLinear()

function Gridap.Arrays.array_cache(a::CellBoundaryVectorFromFacetVector{T}) where T
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

function Gridap.Arrays.getindex!(cache,a::CellBoundaryVectorFromFacetVector,cell::Integer)
  fvc,cwfc,vbc=cache
  vb=vbc[cell]
  cwf=getindex!(cwfc,a.cell_wise_facets_ids,cell)
  for (lfacet,gfacet) in enumerate(cwf)
    fv=getindex!(fvc[cell][lfacet],a.facet_vector,gfacet)
    vb.array[lfacet]=fv
  end
  vb
end

function Base.getindex(a::CellBoundaryVectorFromFacetVector,cell::Integer)
  c=array_cache(a)
  Gridap.Arrays.getindex!(c,a,cell)
end

struct CellBoundaryOpt{M<:DiscreteModel,TBT<:Triangulation,SF}
  model::M
  btrian::TBT
  sign_flip::SF
  function CellBoundaryOpt(model::M) where M<:DiscreteModel
    face_to_bgface = collect(1:num_facets(model))
    btrian=BoundaryTriangulation(model,
                              face_to_bgface,
                              Fill(Int8(1),num_facets(model)))
    TBT=typeof(btrian)
    # TO-DO: here I am reusing the machinery for global RT FE spaces.
    #        Sure there is a way to decouple this from global RT FE spaces.
    function _get_sign_flip(model)
      basis,reffe_args,reffe_kwargs = ReferenceFE(raviart_thomas,Float64,0)
      cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
      Gridap.FESpaces.get_sign_flip(model,cell_reffe)
    end
    sign_flip=_get_sign_flip(model)
    SF=typeof(sign_flip)
    new{M,TBT,SF}(model,btrian,sign_flip)
  end
end


struct CellBoundaryOwnerNref{T,G,SF} <: AbstractVector{Gridap.Fields.VectorBlock{T}}
    cell_boundary_nref::CellBoundaryCompressedVector{T,G}
    sign_flip::SF
end

function Gridap.Arrays.array_cache(a::CellBoundaryOwnerNref{T}) where {T}
  array_cache(a.cell_boundary_nref),array_cache(a.sign_flip)
end

function Base.getindex(a::CellBoundaryOwnerNref,cell::Integer)
 c=array_cache(a)
 Gridap.Arrays.getindex!(c,a,cell)
end

function Gridap.Arrays.getindex!(cache,a::CellBoundaryOwnerNref,cell::Integer)
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

Base.size(a::CellBoundaryOwnerNref) = size(a.cell_boundary_nref)
Base.IndexStyle(::Type{<:CellBoundaryOwnerNref}) = IndexLinear()

function _cell_lface_to_nref(cb::CellBoundaryOpt)
  glue = cb.btrian.glue
  cell_trian = cb.btrian.cell_trian
  ## Reference normal
  function f(r)
    p = Gridap.ReferenceFEs.get_polytope(r)
    lface_to_n = Gridap.ReferenceFEs.get_facet_normal(p)
    lface_to_pindex_to_perm = Gridap.ReferenceFEs.get_face_vertex_permutations(p,num_cell_dims(p)-1)
    nlfaces = length(lface_to_n)
    lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
    lface_pindex_to_n
  end
  ctype_lface_pindex_to_nref = map(f, get_reffes(cell_trian))
  CellBoundaryCompressedVector(ctype_lface_pindex_to_nref,glue)
end

function _cell_lface_to_owner_nref(cb::CellBoundaryOpt)
  cell_lface_to_nref=_cell_lface_to_nref(cb)
  CellBoundaryOwnerNref(cell_lface_to_nref,cb.sign_flip)
end

function get_cell_normal_vector(cb::CellBoundaryOpt)
  _get_cell_normal_vector(cb,_cell_lface_to_nref)
end

function get_cell_owner_normal_vector(cb::CellBoundaryOpt)
  _get_cell_normal_vector(cb,_cell_lface_to_owner_nref)
end

function _get_cell_normal_vector(cb::CellBoundaryOpt,cell_lface_to_nref::Function)
  glue = cb.btrian.glue
  cell_trian = cb.btrian.cell_trian

  cell_lface_to_nref = cell_lface_to_nref(cb)
  cell_lface_s_nref = lazy_map(Gridap.Fields.constant_field,cell_lface_to_nref)

  # Inverse of the Jacobian transpose
  cell_q_x = get_cell_map(cell_trian)
  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(Gridap.Fields.pinvJt),cell_q_Jt)
  cell_lface_q_invJt = transform_cell_to_cell_lface_array(glue, cell_q_invJt)

  # Change of domain
  cell_lface_s_q = get_cell_ref_map(cb)
  cell_lface_s_invJt = lazy_map(∘,cell_lface_q_invJt,cell_lface_s_q)
  #face_s_n =
  lazy_map(Broadcasting(Operation(Gridap.Geometry.push_normal)),
           cell_lface_s_invJt,
           cell_lface_s_nref)
  #Fields.MemoArray(face_s_n)
end

function Gridap.Geometry.get_cell_map(cb::CellBoundaryOpt)
  facet_map=get_cell_map(cb.btrian)
  CellBoundaryVectorFromFacetVector(cb.btrian.glue,
                                    _get_cell_wise_facets(cb),
                                    facet_map)
  # cell_map=get_cell_map(cb.btrian.cell_trian)
  # cell_lface_ref_map=get_cell_ref_map(cb)
  # r=lazy_map(Broadcasting(∘),cell_map,cell_lface_ref_map)
end

function _get_cell_wise_facets(cb::CellBoundaryOpt)
  model = cb.model
  gtopo = get_grid_topology(model)
  D     = num_cell_dims(model)
  Gridap.Geometry.get_faces(gtopo, D, D-1)
end

function Gridap.Geometry.get_cell_ref_map(cb::CellBoundaryOpt)
  cell_lface_to_q_vertex_coords = _compute_cell_lface_to_q_vertex_coords(cb)
  f(p) = Gridap.ReferenceFEs.get_shapefuns(Gridap.ReferenceFEs.LagrangianRefFE(Float64,Gridap.ReferenceFEs.get_polytope(p),1))
  ftype_to_shapefuns = map( f, Gridap.Geometry.get_reffes(cb.btrian) )
  face_to_shapefuns = expand_cell_data(ftype_to_shapefuns,cb.btrian.glue.face_to_ftype)
  cell_to_lface_to_shapefuns = transform_face_to_cell_lface_array(cb.btrian.glue,face_to_shapefuns)
  face_s_q = lazy_map(Gridap.Fields.linear_combination,
                      cell_lface_to_q_vertex_coords,
                      cell_to_lface_to_shapefuns)
end

function _compute_cell_lface_to_q_vertex_coords(cb::CellBoundaryOpt)
  ctype_to_lface_to_pindex_to_qcoords=Gridap.Geometry._compute_face_to_q_vertex_coords(cb.btrian)
  CellBoundaryCompressedVector(
    ctype_to_lface_to_pindex_to_qcoords.ctype_lface_pindex_to_value,
    cb.btrian.glue)
end

function transform_face_to_cell_lface_array(glue,
                                            face_array::Gridap.Arrays.CompressedArray,
                                            f::Function=identity)
  T=typeof(f(face_array.values[1]))
  ctype_to_vector_block=
    Vector{Gridap.Fields.VectorBlock{T}}(undef,length(glue.ctype_to_lface_to_ftype))
  for ctype=1:length(glue.ctype_to_lface_to_ftype)
     num_facets=length(glue.ctype_to_lface_to_ftype[ctype])
     v=Vector{T}(undef,num_facets)
     t=Vector{Bool}(undef,num_facets)
     t.=true
     for lface=1:num_facets
      ftype=glue.ctype_to_lface_to_ftype[ctype][lface]
      v[lface]=f(face_array.values[ftype])
     end
     ctype_to_vector_block[ctype]=Gridap.Fields.ArrayBlock(v,t)
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block,glue.cell_to_ctype)
end

function transform_face_to_cell_lface_expanded_array(
  glue,
  face_array::Gridap.Arrays.LazyArray{<:Fill{typeof(transpose)}})
  Gridap.Helpers.@check typeof(face_array.args[1]) <: Gridap.Arrays.CompressedArray

  T = Gridap.Arrays.return_type(transpose,face_array.args[1].values[1])
  v = Vector{T}(undef,length(face_array.args[1].values))
  for i=1:length(face_array.args[1].values)
    v[i]=evaluate(transpose,face_array.args[1].values[i])
  end
  face_array_compressed=Gridap.Arrays.CompressedArray(v,face_array.args[1].ptrs)
  transform_face_to_cell_lface_expanded_array(glue,face_array_compressed)
end


function transform_face_to_cell_lface_expanded_array(glue,
                                                     face_array::Gridap.Arrays.CompressedArray)
  ftype_to_block_layout=_get_block_layout(face_array.values)
  T=eltype(face_array.values[1])
  if length(ftype_to_block_layout[1][1]) == 1
    TB=Gridap.Fields.VectorBlock{T}
    TF1=Gridap.Fields.VectorBlock{TB}
  else
    Gridap.Helpers.@check length(ftype_to_block_layout[1][1])==2
    TB=Gridap.Fields.MatrixBlock{T}
    TF1=Gridap.Fields.MatrixBlock{TB}
  end
  ctype_to_vector_block= #[c][f1][b][f2] or [c][f1][1,b][1,f2]
    Vector{Gridap.Fields.VectorBlock{TF1}}(undef,length(glue.ctype_to_lface_to_ftype))
  for ctype=1:length(glue.ctype_to_lface_to_ftype)
     num_facets=length(glue.ctype_to_lface_to_ftype[ctype])
     vf1=Vector{TF1}(undef,num_facets)
     tf1=Vector{Bool}(undef,num_facets)
     tf1.=true
     for lface=1:num_facets
      ftype=glue.ctype_to_lface_to_ftype[ctype][lface]
      if length(ftype_to_block_layout[ftype][1])==1
        vb = Vector{TB}(undef,length(face_array.values[ftype]))
        tb = Vector{Bool}(undef,length(face_array.values[ftype]))
      else
        vb = Matrix{TB}(undef,(1,length(face_array.values[ftype])))
        tb = Matrix{Bool}(undef,(1,length(face_array.values[ftype])))
      end
      tb .= false
      for blk=1:length(face_array.values[ftype])
        if face_array.values[ftype].touched[blk]
          if length(ftype_to_block_layout[ftype][1])==1
            vf2 = Vector{T}(undef,num_facets)
            tf2 = Vector{Bool}(undef,num_facets)
            tf2.= false
            tf2[lface]=true
            vf2[lface]=face_array.values[ftype].array[blk]
            vb[blk]=Gridap.Fields.ArrayBlock(vf2,tf2)
            tb[blk]=true
          else
            Gridap.Helpers.@check length(ftype_to_block_layout[ftype][1])==2
            vf2 = Matrix{T}(undef,(1,num_facets))
            tf2 = Matrix{Bool}(undef,(1,num_facets))
            tf2.= false
            tf2[1,lface]=true
            vf2[1,lface]=face_array.values[ftype].array[1,blk]
            vb[1,blk]=Gridap.Fields.ArrayBlock(vf2,tf2)
            tb[1,blk]=true
          end
        end
      end
      vf1[lface]=Gridap.Fields.ArrayBlock(vb,tb)
     end
     ctype_to_vector_block[ctype]=Gridap.Fields.ArrayBlock(vf1,tf1)
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block,glue.cell_to_ctype)
end


function transform_cell_to_cell_lface_array(glue,
                                            cell_array::Fill)
    d = Gridap.Arrays.CompressedArray([cell_array.value,],Fill(1,length(cell_array)))
    transform_cell_to_cell_lface_array(glue,d)
end


function transform_cell_to_cell_lface_array(glue,
                                            cell_array::Gridap.Arrays.CompressedArray)
  T=typeof(cell_array.values[1])
  ctype_to_vector_block=
    Vector{Gridap.Fields.VectorBlock{T}}(undef,length(glue.ctype_to_lface_to_ftype))
  for ctype=1:length(glue.ctype_to_lface_to_ftype)
     num_facets=length(glue.ctype_to_lface_to_ftype[ctype])
     v=Vector{T}(undef,num_facets)
     t=Vector{Bool}(undef,num_facets)
     t.=true
     for lface=1:num_facets
      v[lface]=cell_array.values[ctype]
     end
     ctype_to_vector_block[ctype]=Gridap.Fields.ArrayBlock(v,t)
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block,glue.cell_to_ctype)
end

function transform_cell_to_cell_lface_array(
  glue,
  cell_array::Gridap.Arrays.LazyArray{<:Fill{typeof(transpose)}})
  Gridap.Helpers.@check typeof(cell_array.args[1]) <: Fill
  cell_array_fill = Fill(evaluate(transpose,cell_array.args[1].value),length(cell_array))
  transform_cell_to_cell_lface_array(glue,cell_array_fill)
end

function Gridap.Arrays.lazy_map(k::Gridap.Fields.LinearCombinationMap,
                  ::Type{T},b::CellBoundaryCompressedVector,c::Fill) where T
  d = Gridap.Arrays.CompressedArray([c.value,],Fill(1,length(c)))
  lazy_map(k,T,b,d)
end

function Gridap.Arrays.lazy_map(k::Gridap.Fields.LinearCombinationMap,
                                ::Type{T},
                                b::CellBoundaryCompressedVector,
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
    CellBoundaryCompressedVector(ctype_lface_pindex_to_r,b.glue)
  else
    Gridap.Helpers.@notimplemented
  end
end

function Gridap.Arrays.lazy_map(k::typeof(evaluate),
                                ::Type{T},
                                b::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock},
                                c::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock}) where T
  Gridap.Helpers.@check length(b) == length(c)
  if c.ptrs === b.ptrs || c.ptrs == b.ptrs
     values_r = Vector{T}(undef,length(b.values))
     for i=1:length(b.values)
        values_r[i] =evaluate(k,b.values[i],c.values[i])
     end
  else
    Gridap.Helpers.@notimplemented
  end
  Gridap.Arrays.CompressedArray(values_r,b.ptrs)
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



function Gridap.Arrays.lazy_map(
  k::typeof(Gridap.Arrays.evaluate),
  ::Type{T},
  a::Fill,
  b::CellBoundaryCompressedVector{<:AbstractArray{<:Point}}) where T
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
  CellBoundaryCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

function Gridap.Arrays.lazy_map(
  k::typeof(evaluate),
  ::Type{T},
  a::Gridap.Arrays.CompressedArray,
  b::CellBoundaryCompressedVector) where T

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
  CellBoundaryCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

function Gridap.Arrays.lazy_map(
  k::typeof(evaluate),
  ::Type{T},
  b::CellBoundaryCompressedVector,
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
  CellBoundaryCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

# function lazy_map(
#   k::typeof(evaluate),
#   ::Type{T},
#   a::Fill,
#   b::CellBoundaryCompressedVector,
#   c::CellBoundaryCompressedVector) where T
#   if b.glue !== c.glue
#     return LazyArray(T,a,b,c)
#   end
#   @check length(a) == length(b)
#   ctype_lface_pindex_to_r = Vector{Vector{Vector{T}}}(undef,length(b.ctype_lface_pindex_to_value))
#   for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
#     lface_pindex_to_r = Vector{Vector{T}}(undef,length(lface_pindex_to_value))
#     for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
#       pindex_to_r = Vector{T}(undef,length(pindex_to_value))
#       ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
#       if ftype != UNSET
#         for (pindex, bvalue) in enumerate(pindex_to_value)
#           cvalue = c.ctype_lface_pindex_to_value[ctype][lface][pindex]
#           pindex_to_r[pindex] = evaluate(a.value,bvalue,cvalue)
#         end
#       end
#       lface_pindex_to_r[lface] = pindex_to_r
#     end
#     ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
#   end
#   CellBoundaryCompressedVector(ctype_lface_pindex_to_r,b.glue)
# end

function Gridap.Arrays.lazy_map(
  k::typeof(evaluate),
  ::Type{T},
  b::CellBoundaryVectorFromFacetVector,
  a::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock}) where T
  Gridap.Helpers.@check length(a) == length(b)
  af=_cell_lfacet_vector_to_facet_vector(b.glue,a)
  bf=b.facet_vector
  bfx=lazy_map(evaluate,bf,af)
  CellBoundaryVectorFromFacetVector(b.glue,b.cell_wise_facets_ids,bfx)
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
  b::CellBoundaryVectorFromFacetVector)
  ∇bf=lazy_map(∇,b.facet_vector)
  CellBoundaryVectorFromFacetVector(b.glue,b.cell_wise_facets_ids,∇bf)
end

function quadrature_points_and_weights(cb::CellBoundaryOpt, degree::Integer)
  model = cb.model
  D = num_cell_dims(model)
  p = lazy_map(Gridap.ReferenceFEs.get_polytope,get_reffes(model))
  Gridap.Helpers.@check length(p) == 1

  f(p) = Gridap.ReferenceFEs.Quadrature(Gridap.ReferenceFEs.get_polytope(p),degree)
  ftype_to_quads = map( f, Gridap.Geometry.get_reffes(cb.btrian) )
  face_to_quads = expand_cell_data(ftype_to_quads,cb.btrian.glue.face_to_ftype)

  xf_array_block = transform_face_to_cell_lface_array(
         cb.btrian.glue,face_to_quads,Gridap.ReferenceFEs.get_coordinates)

  wf_array_block = transform_face_to_cell_lface_array(
    cb.btrian.glue,face_to_quads,Gridap.ReferenceFEs.get_weights)

  (xf_array_block,wf_array_block)
end

# Optimization for
#
#  g = lazy_map( linear_combination, cell_to_i_to_val,cell_to_i_to_f)
#  lazy_map(evaluate,g,cell_to_x)
#
# Required for the following lines in RTH driver
#    ∂Topt=CellBoundaryOpt(model)
#    m=Gridap.Geometry.get_cell_ref_map(∂Topt)
#    xopt,wopt=quadrature_points_and_weights(∂Topt,2)
#    mx=lazy_map(evaluate,m,xopt)

function Gridap.Arrays.lazy_map(
  ::typeof(evaluate), a::Gridap.Arrays.LazyArray{<:Fill{typeof(Gridap.Fields.linear_combination)}},
                      x::Gridap.Arrays.CompressedArray{<:Gridap.Fields.VectorBlock})

  i_to_values = a.args[1]
  i_to_basis = a.args[2]
  i_to_basis_x = lazy_map(evaluate,i_to_basis,x)
  lazy_map(Gridap.Fields.LinearCombinationMap(:),i_to_values,i_to_basis_x)
end

function Gridap.Arrays.return_cache(
  #::typeof(evaluate),
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
  #::typeof(evaluate),
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


function Gridap.Arrays.return_cache(
  a::Gridap.Fields.Field,
  x::Gridap.Fields.VectorBlock)

   Gridap.Helpers.@check all(x.touched)
   T=Gridap.Arrays.return_type(evaluate,a,x.array[1])
   Tc=typeof(Gridap.Arrays.return_cache(evaluate,a,x.array[1]))
   r=Vector{T}(undef,length(x.array))
   rc=Vector{Tc}(undef,length(x.array))
   rc[1]=Gridap.Arrays.return_cache(evaluate,a,x.array[1])
   for i=2:length(x.array)
      if x.array[i] === x.array[1]
        rc[i]=rc[1]
      else
        rc[i]=Gridap.Arrays.return_cache(evaluate,a,x.array[i])
      end
   end
   (Gridap.Fields.ArrayBlock(r,x.touched),Gridap.Fields.ArrayBlock(rc,x.touched))
end

function Gridap.Arrays.evaluate!(
  cache,
  a::Gridap.Fields.Field,
  x::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check all(x.touched)
  (r,rc)=cache
  for i=1:length(x)
    r.array[i]=Gridap.Arrays.evaluate!(rc.array[i],evaluate,a,x.array[i])
  end
  r
end

function restrict_to_cell_boundary(cb::CellBoundaryOpt,
                                   fe_basis::Union{Gridap.FESpaces.FEBasis,Gridap.CellData.CellField})
  model = cb.model
  D = num_cell_dims(model)
  trian = get_triangulation(fe_basis)
  if isa(trian,Triangulation{D-1,D})
    _restrict_to_cell_boundary_facet_fe_basis(cb,fe_basis)
  elseif isa(trian,Triangulation{D,D})
    _restrict_to_cell_boundary_cell_fe_basis(cb,fe_basis)
  end
end

function _restrict_to_cell_boundary_cell_fe_basis(cb::CellBoundaryOpt,
                                                  cell_fe_basis::Union{Gridap.FESpaces.FEBasis,Gridap.CellData.CellField})
  # TO-THINK:
  #     1) How to deal with CellField objects which are NOT FEBasis objects?
  #     2) How to deal with differential operators applied to FEBasis objects?
  D = num_cell_dims(cb.model)
  Gridap.Helpers.@check isa(get_triangulation(cell_fe_basis),Triangulation{D,D})
  cell_a_q = transform_cell_to_cell_lface_array(cb.btrian.glue,Gridap.CellData.get_data(cell_fe_basis))
  cell_s2q = get_cell_ref_map(cb)
  lazy_map(Broadcasting(∘),cell_a_q,cell_s2q)
end

function _restrict_to_cell_boundary_facet_fe_basis(cb::CellBoundaryOpt,
                                                   facet_fe_basis::Gridap.FESpaces.FEBasis)

  # TO-THINK:
  #     1) How to deal with CellField objects which are NOT FEBasis objects?
  #     2) How to deal with differential operators applied to FEBasis/CellField objects?
  D = num_cell_dims(cb.model)
  Gridap.Helpers.@check isa(get_triangulation(facet_fe_basis),Triangulation{D-1,D})

  transform_face_to_cell_lface_expanded_array(
    cb.btrian.glue,
    Gridap.CellData.get_data(facet_fe_basis))
end


function Gridap.Arrays.return_value(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.VectorBlock,
  h::Gridap.Fields.VectorBlock)
  evaluate(k,f,h)
end

function Gridap.Arrays.return_cache(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.VectorBlock,
  h::Gridap.Fields.VectorBlock) where {A,N}
  Gridap.Helpers.@check length(f)==length(h)
  Gridap.Helpers.@check f.touched==h.touched
  fi = Gridap.Arrays.testitem(f)
  hi = Gridap.Arrays.testitem(h)
  li = Gridap.Arrays.return_cache(k,fi,hi)
  fix = Gridap.Arrays.evaluate!(li,k,fi,hi)
  l = Vector{typeof(li)}(undef,size(h.array))
  g = Vector{typeof(fix)}(undef,size(h.array))
  for i in eachindex(h.array)
    if h.touched[i]
      l[i] = Gridap.Arrays.return_cache(k,f.array[i],h.array[i])
    end
  end
  Gridap.Fields.ArrayBlock(g,h.touched),l
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.VectorBlock,
  h::Gridap.Fields.VectorBlock)

  g,l = cache
  Gridap.Helpers.@check g.touched == h.touched
  for i in eachindex(h.array)
    if h.touched[i]
      g.array[i] = Gridap.Arrays.evaluate!(l[i],k,f.array[i],h.array[i])
    end
  end
  g
end

function Gridap.Arrays.return_value(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.Field,
  h::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}) where {N}
  evaluate(k,f,h)
end

function Gridap.Arrays.return_cache(
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.Field,
  h::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}) where {N}
  hi = Gridap.Arrays.testitem(h)
  li = Gridap.Arrays.return_cache(k,f,hi)
  fix = Gridap.Arrays.evaluate!(li,k,f,hi)
  l = Array{typeof(li),N}(undef,size(h.array))
  g = Array{typeof(fix),N}(undef,size(h.array))
  for i in eachindex(h.array)
    if h.touched[i]
      l[i] = Gridap.Arrays.return_cache(k,f,h.array[i])
    end
  end
  Gridap.Fields.ArrayBlock(g,h.touched),l
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Broadcasting{typeof(∘)},
  f::Gridap.Fields.Field,
  h::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}) where {N}
  g,l = cache
  Gridap.Helpers.@check g.touched == h.touched
  for i in eachindex(h.array)
    if h.touched[i]
      g.array[i] = Gridap.Arrays.evaluate!(l[i],k,f,h.array[i])
    end
  end
  g
end

function Gridap.Fields.constant_field(a::Gridap.Fields.ArrayBlock{T,N}) where {T<:Number,N}
  v=Array{Gridap.Fields.ConstantField{T},N}(undef,size(a))
  for (i,e) in enumerate(a.array)
    v[i]=Gridap.Fields.ConstantField(e)
  end
  Gridap.Fields.ArrayBlock(v,a.touched)
end

function Gridap.Arrays.return_value(
  op::Broadcasting{<:Operation},
  x::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}...) where N
  xi=map(a->a.array[1],x)
  T=Gridap.Arrays.return_type(op,xi...)
  v=Array{T,N}(undef,size(x[1]))
  for i in eachindex(v)
    xi=map(a->a.array[i],x)
    v[i]=Gridap.Arrays.return_value(op,xi...)
  end
  Gridap.Fields.ArrayBlock(v,x[1].touched)
end

function Gridap.Arrays.evaluate!(cache,
  op::Broadcasting{<:Operation},
  x::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}...) where N
  Gridap.Arrays.return_value(op,x...)
end

function Gridap.Arrays.return_cache(
  op::Broadcasting{<:Operation},
  x::Gridap.Fields.ArrayBlock{<:Gridap.Fields.Field,N}...) where N
  xi=map(a->a.array[1],x)
  T=Gridap.Arrays.return_type(op,xi...)
  Tc=typeof(Gridap.Arrays.return_cache(op,xi...))
  v=Array{T,N}(undef,size(x[1]))
  vc=Array{Tc,N}(undef,size(x[1]))
  for i in eachindex(v)
    xi=map(a->a.array[i],x)
    vc[i]=Gridap.Arrays.return_cache(op,xi...)
  end
  Gridap.Fields.ArrayBlock(v,x[1].touched), Gridap.Fields.ArrayBlock(vc,x[1].touched)
end

function Base.:(∘)(a::Gridap.Fields.VectorBlock{<:Gridap.Fields.Field},
                b::Gridap.Fields.VectorBlock{<:Gridap.Fields.Field})
  Gridap.Helpers.@check size(a)==size(b)
  Gridap.Helpers.@check all(a.touched==b.touched)
  Gridap.Helpers.@check all(a.touched)
  T=Gridap.Arrays.return_type(∘,a.array[1],b.array[1])
  v=Vector{T}(undef,length(a.array))
  for i=1:length(v)
    v[i]=a.array[i]∘b.array[i]
  end
  Gridap.Fields.ArrayBlock(v,a.touched)
end

function Gridap.Arrays.return_cache(
  k::Gridap.Fields.IntegrationMap,
  aq::Gridap.Fields.VectorBlock{A},
  w::Gridap.Fields.VectorBlock{B},
  jq::Gridap.Fields.VectorBlock{C}) where{A,B,C}
  Gridap.Helpers.@check length(aq)==length(w)
  Gridap.Helpers.@check length(w)==length(jq)
  Gridap.Helpers.@check all(aq.touched)
  Gridap.Helpers.@check all(w.touched)
  Gridap.Helpers.@check all(jq.touched)
  aqi = Gridap.Arrays.testvalue(A)
  wi  = Gridap.Arrays.testvalue(B)
  jqi = Gridap.Arrays.testvalue(C)
  ci = Gridap.Arrays.return_cache(k,aqi,wi,jqi)
  hi = Gridap.Arrays.evaluate!(ci,k,aqi,wi,jqi)
  a = Vector{typeof(hi)}(undef,size(aq.array))
  b = Vector{typeof(ci)}(undef,size(aq.array))
  for i in eachindex(aq.array)
    b[i] = Gridap.Arrays.return_cache(k,
                        aq.array[i],
                        w.array[i],
                        jq.array[i])
  end
  Gridap.Fields.ArrayBlock(a,aq.touched), b
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Gridap.Fields.IntegrationMap,
  aq::Gridap.Fields.VectorBlock,
  w::Gridap.Fields.VectorBlock,
  jq::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check length(aq)==length(w)
  Gridap.Helpers.@check length(w)==length(jq)
  Gridap.Helpers.@check all(aq.touched)
  Gridap.Helpers.@check all(w.touched)
  Gridap.Helpers.@check all(jq.touched)
  a,b = cache
  # TO-DO: variable number of faces per cell boundary
  Gridap.Helpers.@check length(a)==length(aq)
  for i in eachindex(aq.array)
    a.array[i] = evaluate!(b[i],
                           k,
                           aq.array[i],
                           w.array[i],
                           jq.array[i])
  end
  a
end

struct SumFacetsMap <: Gridap.Fields.Map end

function _sum_facets(cb::CellBoundaryOpt,vx,w,jx)
   int=lazy_map(Gridap.Fields.IntegrationMap(),vx,w,jx)
   lazy_map(SumFacetsMap(),int)
end

function Gridap.Arrays.return_cache(
  k::SumFacetsMap,
  a::Gridap.Fields.VectorBlock)
  Gridap.Helpers.@check all(a.touched)
  m=Gridap.Fields.BroadcastingFieldOpMap(+)
  c=Gridap.Arrays.return_cache(m,a.array[1],a.array[2])
  v=Vector{typeof(c)}(undef,length(a.array)-1)
  v[1]=c
  res=Gridap.Arrays.evaluate!(v[1],m,a.array[1],a.array[2])
  for i=3:length(a.array)
    v[i-1]=Gridap.Arrays.return_cache(m,res,a.array[i])
    res=Gridap.Arrays.evaluate!(v[i-1],m,res,a.array[i])
  end
  v
end

function Gridap.Arrays.evaluate!(
  cache,
  k::SumFacetsMap,
  a::Gridap.Fields.VectorBlock{A}) where{A}
  m=Gridap.Fields.BroadcastingFieldOpMap(+)
  Gridap.Helpers.@check all(a.touched)
  res=Gridap.Arrays.evaluate!(cache[1],m,a.array[1],a.array[2])
  for i=3:length(a.array)
    res=Gridap.Arrays.evaluate!(cache[i-1],m,res,a.array[i])
  end
  res
end


struct RestrictFacetDoFsToCellBoundary{F} <: Gridap.Fields.Map
  facet_dofs::F
end

function Gridap.Arrays.return_cache(
  k::RestrictFacetDoFsToCellBoundary,
  cell_facets::AbstractArray{<:Integer})
  array_cache(k.facet_dofs), Gridap.Arrays.CachedVector(eltype(cell_facets))
end

function Gridap.Arrays.evaluate!(
  cache,
  k::RestrictFacetDoFsToCellBoundary,
  cell_facets::AbstractArray{<:Integer})
  fdc,c=cache
  Gridap.Arrays.setsize!(c,(_count_dofs_cell(cell_facets,k.facet_dofs,fdc),))
  _fill_dofs_cell!(c.array,cell_facets,k.facet_dofs,fdc)
  c.array
end

function _count_dofs_cell(cell_facets,facet_dofs,facet_dofs_cache)
  count=0
  for facet_id in cell_facets
      current_facet_dofs=getindex!(facet_dofs_cache,facet_dofs,facet_id)
      count=count+length(current_facet_dofs)
  end
  count
end

function _fill_dofs_cell!(cell_dofs,cell_facets,facet_dofs,facet_dofs_cache)
  current=1
  for facet_id in cell_facets
      current_facet_dofs=getindex!(facet_dofs_cache,facet_dofs,facet_id)
      for i in current_facet_dofs
        cell_dofs[current]=i
        current=current+1
      end
  end
  current
end

function restrict_facet_dof_ids_to_cell_boundary(cb::CellBoundaryOpt,facet_dof_ids)
  cell_wise_facets = _get_cell_wise_facets(cb)
  m=RestrictFacetDoFsToCellBoundary(facet_dof_ids)
  lazy_map(m,cell_wise_facets)
end
