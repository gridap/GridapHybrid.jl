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
  glue::G
  cell_wise_facets_ids::C
  facet_vector::V
  function CellBoundaryVectorFromFacetVector(glue,cell_wise_facets_ids,facet_vector)
    G=typeof(glue)
    C=typeof(cell_wise_facets_ids)
    V=typeof(facet_vector)
    T=eltype(V)
    new{T,G,C,V}(glue,cell_wise_facets_ids,facet_vector)
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

struct CellBoundary{M<:DiscreteModel,TBT<:Triangulation,SF}
  model::M
  btrian::TBT
  sign_flip::SF
  function CellBoundary(model::M) where M<:DiscreteModel
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

function _cell_lface_to_nref(args...)
  model,glue = args[1],first(args[2:end])
  cell_grid = get_grid(model)
  ## Reference normal
  function f(r)
    p = Gridap.ReferenceFEs.get_polytope(r)
    lface_to_n = Gridap.ReferenceFEs.get_facet_normal(p)
    lface_to_pindex_to_perm = Gridap.ReferenceFEs.get_face_vertex_permutations(p,num_cell_dims(p)-1)
    nlfaces = length(lface_to_n)
    lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
    lface_pindex_to_n
  end
  ctype_lface_pindex_to_nref = map(f, get_reffes(cell_grid))
  CellBoundaryCompressedVector(ctype_lface_pindex_to_nref,glue)
end

function _cell_lface_to_owner_nref(args...)
  model,glue,sign_flip = args
  cell_lface_to_nref=_cell_lface_to_nref(model,glue)
  CellBoundaryOwnerNref(cell_lface_to_nref,sign_flip)
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell.
"""
function get_cell_normal_vector(cb::CellBoundary)
  _get_cell_normal_vector(cb.model, cb.btrian.glue, _cell_lface_to_nref)
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell owner of the facet.
"""
function get_cell_owner_normal_vector(cb::CellBoundary)
  _get_cell_normal_vector(cb.model, cb.btrian.glue, _cell_lface_to_owner_nref, cb.sign_flip)
end

function _get_cell_normal_vector(model,glue,cell_lface_to_nref::Function,sign_flip=nothing)
  cell_grid = get_grid(model)

  cell_lface_to_nref = cell_lface_to_nref(model,glue,sign_flip)
  cell_lface_s_nref = lazy_map(Gridap.Fields.constant_field,cell_lface_to_nref)

  # Inverse of the Jacobian transpose
  cell_q_x = get_cell_map(cell_grid)
  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(Gridap.Fields.pinvJt),cell_q_Jt)
  cell_lface_q_invJt = transform_cell_to_cell_lface_array(glue, cell_q_invJt)

  # Change of domain
  cell_lface_s_q = _setup_tcell_lface_mface_map(num_cell_dims(model)-1,model,glue)

  cell_lface_s_invJt = lazy_map(∘,cell_lface_q_invJt,cell_lface_s_q)
  #face_s_n =
  lazy_map(Broadcasting(Operation(Gridap.Geometry.push_normal)),
           cell_lface_s_invJt,
           cell_lface_s_nref)
  #Fields.MemoArray(face_s_n)
end

function _setup_tcell_lface_mface_map(d,model,glue)
  ctype_to_lface_to_pindex_to_qcoords=Gridap.Geometry._compute_face_to_q_vertex_coords_body(d,model,glue)
  cell_lface_to_q_vertex_coords = CellBoundaryCompressedVector(
                                  ctype_to_lface_to_pindex_to_qcoords.ctype_lface_pindex_to_value,
                                  glue)
  f(p) = Gridap.ReferenceFEs.get_shapefuns(
         Gridap.ReferenceFEs.LagrangianRefFE(Float64,Gridap.ReferenceFEs.get_polytope(p),1))

  ################ TO-IMPROVE
  cell_grid = get_grid(model)
  D = num_cell_dims(model)
  ctype_reffe = get_reffes(cell_grid)
  Gridap.Helpers.@notimplementedif length(ctype_reffe) != 1
  reffe = first(ctype_reffe)
  freffes = get_reffaces(ReferenceFE{D-1},reffe)
  Gridap.Helpers.@notimplementedif length(freffes) != 1
  ftrian_reffes= Fill(first(freffes),length(glue.face_to_ftype))
  ################ TO-IMPROVE

  ftype_to_shapefuns = map( f,  ftrian_reffes)
  face_to_shapefuns = expand_cell_data(ftype_to_shapefuns,glue.face_to_ftype)
  cell_to_lface_to_shapefuns = transform_face_to_cell_lface_array(glue,face_to_shapefuns)
  lazy_map(Gridap.Fields.linear_combination,
           cell_lface_to_q_vertex_coords,
           cell_to_lface_to_shapefuns)
end

function Gridap.Geometry.get_cell_map(cb::CellBoundary)
  facet_map=get_cell_map(cb.btrian)
  CellBoundaryVectorFromFacetVector(cb.btrian.glue,
                                    _get_cell_wise_facets(cb),
                                    facet_map)
  # cell_map=get_cell_map(cb.btrian.cell_trian)
  # cell_lface_ref_map=get_cell_ref_map(cb)
  # r=lazy_map(Broadcasting(∘),cell_map,cell_lface_ref_map)
end

function _get_cell_wise_facets(cb::CellBoundary)
  model = cb.model
  _get_cell_wise_facets(model)
end

function _get_cell_wise_facets(model::DiscreteModel)
  gtopo = get_grid_topology(model)
  D     = num_cell_dims(model)
  Gridap.Geometry.get_faces(gtopo, D, D-1)
end

function _get_cells_around_facets(cb)
  model = cb.model
  _get_cells_around_facets(model)
end

function _get_cells_around_facets(model::DiscreteModel)
  gtopo = get_grid_topology(model)
  D     = num_cell_dims(model)
  Gridap.Geometry.get_faces(gtopo, D-1, D)
end

function _compute_cell_lface_to_q_vertex_coords(cb::CellBoundary)
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
                                                     face_array::Fill)
  T=typeof(face_array.value)
  ctype_to_vector_block=
    Vector{Gridap.Fields.VectorBlock{T}}(undef,length(glue.ctype_to_lface_to_ftype))
  for ctype=1:length(glue.ctype_to_lface_to_ftype)
     num_facets=length(glue.ctype_to_lface_to_ftype[ctype])
     v=Vector{T}(undef,num_facets)
     t=Vector{Bool}(undef,num_facets)
     t.=true
     for lface=1:num_facets
      ftype=glue.ctype_to_lface_to_ftype[ctype][lface]
      v[lface]=face_array.value
     end
     ctype_to_vector_block[ctype]=Gridap.Fields.ArrayBlock(v,t)
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block,glue.cell_to_ctype)
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
                                            cell_array::Fill;
                                            add_naive_innermost_block_level=false)
    d = Gridap.Arrays.CompressedArray([cell_array.value,],Fill(1,length(cell_array)))
    transform_cell_to_cell_lface_array(glue,d;add_naive_innermost_block_level=add_naive_innermost_block_level)
end


function transform_cell_to_cell_lface_array(glue,
                                            cell_array::Gridap.Arrays.CompressedArray;
                                            add_naive_innermost_block_level=false)
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
  if add_naive_innermost_block_level
    ctype_to_vector_block=collect(lazy_map(AddNaiveInnerMostBlockLevelMap(),ctype_to_vector_block))
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block,glue.cell_to_ctype)
end

function transform_cell_to_cell_lface_array(
  glue,
  cell_array::Gridap.Arrays.LazyArray{<:Fill{typeof(transpose)}};
  add_naive_innermost_block_level=false)
  Gridap.Helpers.@check typeof(cell_array.args[1]) <: Fill
  cell_array_fill = Fill(evaluate(transpose,cell_array.args[1].value),length(cell_array))
  transform_cell_to_cell_lface_array(glue,cell_array_fill; add_naive_innermost_block_level=add_naive_innermost_block_level)
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

function quadrature_points_and_weights(cb::CellBoundary, degree::Integer)
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
#    ∂T=CellBoundary(model)
#    m=Gridap.Geometry.get_cell_ref_map(∂T)
#    x,w=quadrature_points_and_weights(∂T,2)
#    mx=lazy_map(evaluate,m,x)

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

function restrict_to_cell_boundary(cb::CellBoundary,
                                   fe_basis::Gridap.CellData.CellField)
  model = cb.model
  D = num_cell_dims(model)
  trian = get_triangulation(fe_basis)
  if isa(trian,Triangulation{D,D})
    tface_to_mface_map = _setup_tcell_lface_mface_map(D-1,model,cb.btrian.glue)
    _restrict_to_cell_boundary_cell_fe_basis(model,
                                              cb.btrian.glue,
                                              tface_to_mface_map,
                                              fe_basis)
  elseif isa(trian,Triangulation{D-1,D})
    _restrict_to_cell_boundary_facet_fe_basis(model,
                                             cb.btrian.glue,
                                             fe_basis)
  end
end

function _restrict_to_cell_boundary_cell_fe_basis(model,
                                                  glue,
                                                  tface_to_mface_map,
                                                  cell_fe_basis::Gridap.CellData.CellField)
  D = num_cell_dims(model)
  Gridap.Helpers.@check isa(get_triangulation(cell_fe_basis),Triangulation{D,D})
  cell_a_q = transform_cell_to_cell_lface_array(glue,
         Gridap.CellData.get_data(cell_fe_basis);
         add_naive_innermost_block_level=true)
  lazy_map(Broadcasting(∘),cell_a_q,tface_to_mface_map)
end

function _restrict_to_cell_boundary_facet_fe_basis(model,
                                                   glue,
                                                   facet_fe_basis::Gridap.CellData.CellField)

  D = num_cell_dims(model)
  Gridap.Helpers.@check isa(get_triangulation(facet_fe_basis),Triangulation{D-1,D})

  transform_face_to_cell_lface_expanded_array(
    glue,
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

struct RestrictFacetDoFsToCellBoundary{F} <: Gridap.Fields.Map
  facet_dofs::F
end

function Gridap.Arrays.return_cache(
  k::RestrictFacetDoFsToCellBoundary,
  cell_facets::AbstractArray{<:Integer})
  if (isa(k.facet_dofs,AbstractVector{<:Number}))
    T=eltype(k.facet_dofs)
  elseif (isa(k.facet_dofs,AbstractVector{<:AbstractVector{<:Number}}))
    T=eltype(eltype(k.facet_dofs))
  else
    @assert false
  end
  array_cache(k.facet_dofs), Gridap.Arrays.CachedVector(T)
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

function restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets,facet_dof_ids)
  m=RestrictFacetDoFsToCellBoundary(facet_dof_ids)
  lazy_map(m,cell_wise_facets)
end

function _get_block_layout(fields_array::AbstractArray{<:AbstractArray{<:Gridap.Fields.Field}})
  Fill((1,1),length(fields_array))
end

function _get_block_layout(fields_array::AbstractArray{<:Gridap.Fields.ArrayBlock})
  lazy_map(x->((size(x),findall(x.touched))),fields_array)
end

# This function will eventually play the role of its counterpart in Gridap
# function integrate(f::CellField,quad::CellQuadrature)
# TO-THINK: mh,uh currently as LazyArray{...}
#    Not sure if we are loosing generality by constraining a to be of type
#    LazyArray{...}. I need to restrict the cell-wise block array to each
#    individual block, and with a LazyArray{...} this is very efficient as
#    the array is already restricted to each block in the a.args member variable.

#∫( mh*(uh⋅n) )*dK
function integrate_mh_mult_uh_cdot_n_low_level(cb,
  mh,
  uh,
  x::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractArray{<:Point}}},
  w::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractVector}})

  n=get_cell_normal_vector(cb)

  # (uh⋅n)
  uhx=lazy_map(evaluate,uh,x)
  nx=lazy_map(evaluate,n,x)
  uhx_cdot_nx=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(⋅), uhx, nx)

  # mh*(uh⋅n)
  mhx=lazy_map(evaluate,mh,x)
  mhx_mult_uhx_cdot_nx = lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), mhx, uhx_cdot_nx )

  j=lazy_map(∇,get_cell_map(cb))
  jx=lazy_map(evaluate,j,x)

  lazy_map(IntegrationMap(),mhx_mult_uhx_cdot_nx,w,jx)
end


#∫( (vh⋅n)*lh )*dK
function integrate_vh_cdot_n_mult_lh_low_level(
  cb,
  vh,
  lh,
  x::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractArray{<:Point}}},
  w::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractVector}})

  n=get_cell_normal_vector(cb)

  # (vh⋅n)
  vhx=lazy_map(evaluate,vh,x)
  nx=lazy_map(evaluate,n,x)
  vhx_cdot_nx=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(⋅), vhx, nx)

  # (vh⋅n)*lh
  lhx=lazy_map(evaluate,lh,x)
  vhx_cdot_nx_mult_lhx = lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), vhx_cdot_nx, lhx)

  j=lazy_map(∇,get_cell_map(cb))
  jx=lazy_map(evaluate,j,x)

  lazy_map(IntegrationMap(),vhx_cdot_nx_mult_lhx,w,jx)
end

# ∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
function integrate_qh_mult_uh_cdot_n_plus_stab_low_level(
  cb,
  qh,
  uh,
  τ ,
  ph,
  lh,
  x::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractArray{<:Point}}},
  w::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractVector}})

  # (n⋅no)
  n=get_cell_normal_vector(cb)
  no=get_cell_owner_normal_vector(cb)
  nx=lazy_map(evaluate,n,x)
  nox=lazy_map(evaluate,no,x)
  nx_cdot_nox=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(⋅), nx, nox)

  # ∫(qh*(uh⋅n)d∂K
  qhx=lazy_map(evaluate,qh,x)
  uhx=lazy_map(evaluate,uh,x)
  uhx_cdot_nx=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(⋅), uhx, nx)
  qhx_mult_uhx_cdot_nx=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), qhx, uhx_cdot_nx)

  # ∫(qh*(τ*(ph)*(n⋅no))d∂K
  phx=lazy_map(evaluate,ph,x)
  phx_mult_n_cdot_no=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), phx, nx_cdot_nox)
  τx=lazy_map(evaluate,τ,x)
  τx_mult_phx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), τx, phx_mult_n_cdot_no)
  qhx_mult_τx_mult_phx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), qhx, τx_mult_phx_mult_n_cdot_no)

  # -∫(qh*(τ*(lh)*(n⋅no))d∂K
  lhx=lazy_map(evaluate,lh,x)
  lhx_mult_n_cdot_no=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), lhx, nx_cdot_nox)
  τx_mult_lhx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), τx, lhx_mult_n_cdot_no)
  qhx_mult_τx_mult_lhx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), qhx, τx_mult_lhx_mult_n_cdot_no)

  j=lazy_map(∇,get_cell_map(cb))
  jx=lazy_map(evaluate,j,x)

  arg1=lazy_map(Broadcasting(+),
                lazy_map(IntegrationMap(),qhx_mult_uhx_cdot_nx,w,jx),
                lazy_map(IntegrationMap(),qhx_mult_τx_mult_phx_mult_n_cdot_no,w,jx))

  arg2=lazy_map(IntegrationMap(),qhx_mult_τx_mult_lhx_mult_n_cdot_no,w,jx)

  lazy_map(Broadcasting(-),arg1,arg2)
end

# ∫(mh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
function integrate_mh_mult_uh_cdot_n_plus_stab_low_level(
  cb,
  mh,
  uh,
  τ,
  ph,
  lh,
  x::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractArray{<:Point}}},
  w::AbstractArray{<:Gridap.Fields.ArrayBlock{<:AbstractVector}})

  # (n⋅no)
  n=get_cell_normal_vector(cb)
  no=get_cell_owner_normal_vector(cb)
  nx=lazy_map(evaluate,n,x)
  nox=lazy_map(evaluate,no,x)
  nx_cdot_nox=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(⋅), nx, nox)

  # ∫(mh*(uh⋅n)d∂K
  mhx=lazy_map(evaluate,mh,x)
  uhx=lazy_map(evaluate,uh,x)
  uhx_cdot_nx=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(⋅), uhx, nx)
  mhx_mult_uhx_cdot_nx=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), mhx, uhx_cdot_nx)

  # ∫(mh*(τ*(ph)*(n⋅no))d∂K
  phx=lazy_map(evaluate,ph,x)
  phx_mult_n_cdot_no=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), phx, nx_cdot_nox)
  τx=lazy_map(evaluate,τ,x)
  τx_mult_phx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), τx, phx_mult_n_cdot_no)
  mhx_mult_τx_mult_phx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), mhx, τx_mult_phx_mult_n_cdot_no)

  # -∫(mh*(τ*(lh)*(n⋅no))d∂K
  mhx=lazy_map(evaluate,mh,x)
  lhx=lazy_map(evaluate,lh,x)
  lhx_mult_n_cdot_no=lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), lhx, nx_cdot_nox)
  τx_mult_lhx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), τx, lhx_mult_n_cdot_no)
  mhx_mult_τx_mult_lhx_mult_n_cdot_no=
      lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*), mhx, τx_mult_lhx_mult_n_cdot_no)

  j=lazy_map(∇,get_cell_map(cb))
  jx=lazy_map(evaluate,j,x)

  arg1=lazy_map(IntegrationMap(),mhx_mult_uhx_cdot_nx,w,jx)
  arg2=lazy_map(IntegrationMap(),mhx_mult_τx_mult_phx_mult_n_cdot_no,w,jx)
  arg3=lazy_map(IntegrationMap(),mhx_mult_τx_mult_lhx_mult_n_cdot_no,w,jx)

  result=lazy_map(Broadcasting(+),arg1,arg2)
  result=lazy_map(Broadcasting(-),result,arg3)
end

function _generate_glue_among_facet_and_cell_wise_dofs_arrays(
  cells_around_facets,
  cell_wise_facets,
  facet_dof_ids)

  c1=array_cache(cells_around_facets)
  c2=array_cache(cell_wise_facets)
  c3=array_cache(facet_dof_ids)

  result=Vector{NTuple{3,Int}}(undef,length(facet_dof_ids))
  current=1
  ndofs=0
  for facet_gid=1:length(cells_around_facets)
    cell_gid=Gridap.Arrays.getindex!(c1,cells_around_facets,facet_gid)[1]
    current_cell_facets=Gridap.Arrays.getindex!(c2,cell_wise_facets,cell_gid)
    pos=1
    for facet_gid_in_cell in current_cell_facets
      ndofs=length(Gridap.Arrays.getindex!(c3,facet_dof_ids,facet_gid_in_cell))
      if (facet_gid == facet_gid_in_cell)
        break
      else
        pos=pos+ndofs
      end
    end
    result[facet_gid]=(cell_gid,pos,ndofs)
  end
  result
end

struct ExtractFacetDofsFromCellDofs{T<:AbstractVector{<:AbstractVector}} <: Gridap.Fields.Map
   cell_dofs::T
end

function Gridap.Arrays.return_cache(k::ExtractFacetDofsFromCellDofs,
                                    cellgid_facetlpos_ndofs::NTuple{3})
  cell_dofs_cache  = Gridap.Arrays.array_cache(k.cell_dofs)
  T=eltype(eltype(k.cell_dofs))
  facet_dofs_cache = Gridap.Arrays.CachedArray(zeros(T,cellgid_facetlpos_ndofs[3]))
  cell_dofs_cache, facet_dofs_cache
end

function Gridap.Arrays.evaluate!(cache,
                                 k::ExtractFacetDofsFromCellDofs,
                                 cellgid_facetlpos_ndofs::NTuple{3})
  cell_dofs_cache, facet_dofs_cache = cache
  cellgid,facetlpos,ndofs=cellgid_facetlpos_ndofs
  Gridap.Arrays.setsize!(facet_dofs_cache,(ndofs,))
  facet_dofs  = facet_dofs_cache.array
  cell_dofs   = Gridap.Arrays.getindex!(cell_dofs_cache,k.cell_dofs,cellgid)
  facet_dofs .= cell_dofs[facetlpos:facetlpos+ndofs-1]
end


function convert_cell_wise_dofs_array_to_facet_dofs_array(
       cells_around_facets,
       cell_wise_facets,
       cell_dofs_array::AbstractVector{<:AbstractVector},
       facet_dof_ids)
  glue = _generate_glue_among_facet_and_cell_wise_dofs_arrays(
    cells_around_facets, cell_wise_facets, facet_dof_ids)
  k=ExtractFacetDofsFromCellDofs(cell_dofs_array)
  lazy_map(k,glue)
end
