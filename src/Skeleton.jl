using Gridap.Helpers
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Arrays
using Gridap.Visualization
using FillArrays

export Skeleton

struct SkeletonGrid{Dc,Dp,P,A,B,C,D} <: Grid{Dc,Dp}
  parent::P
  node_coord::A
  cell_lface_nodes::B
  ftype_freffe::C
  cell_lface_ftype::D
  function SkeletonGrid(grid::Grid)
    D = num_cell_dims(grid)
    Dp = num_point_dims(grid)
    node_coord = get_node_coordinates(grid)
    cell_ctype = get_cell_type(grid)
    ctype_reffe       = get_reffes(grid)
    @notimplementedif length(ctype_reffe) != 1
    reffe = first(ctype_reffe)
    freffes = get_reffaces(ReferenceFE{D-1},reffe)
    @notimplementedif length(freffes) != 1
    ftype_freffe = [first(freffes),]

    ####
    ctype_lface_ftype = map(reffe->get_face_type(reffe,D-1),ctype_reffe)
    cell_lface_ftype  = expand_cell_data(ctype_lface_ftype,cell_ctype)
    ###

    cell_nodes = get_cell_node_ids(grid)
    function f(reffe)
      lface_to_lnodes=get_face_nodes(reffe,D-1)
      ArrayBlock(lface_to_lnodes,[true for i=1:length(lface_to_lnodes)])
    end
    ctype_lface_to_lnodes = map(f,ctype_reffe)
    cell_lface_to_lnodes  = expand_cell_data(ctype_lface_to_lnodes,cell_ctype)
    m=lazy_map(Reindex,cell_nodes)
    fi = testitem(cell_lface_to_lnodes)
    mi = testitem(m)
    T = return_type(mi, fi)
    cell_lface_nodes=LazyArray(T,m,cell_lface_to_lnodes)

    A = typeof(node_coord)
    B = typeof(cell_lface_nodes)
    C = typeof(ftype_freffe)
    E = typeof(cell_lface_ftype)
    P = typeof(grid)
    new{D-1,Dp,P,A,B,C,E}(
      grid,node_coord,cell_lface_nodes,ftype_freffe,cell_lface_ftype)
  end
end

Geometry.get_node_coordinates(a::SkeletonGrid) = a.node_coord
Geometry.get_cell_node_ids(a::SkeletonGrid) = a.cell_lface_nodes
Geometry.get_reffes(a::SkeletonGrid) = a.ftype_freffe
Geometry.get_cell_type(a::SkeletonGrid) = a.cell_lface_ftype


function Geometry.get_cell_coordinates(trian::SkeletonGrid)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_node_ids(trian)
  lazy_map(Reindex(node_to_coords),cell_to_nodes)
end

# The following overloads of Gridap API with Reindex are needed to support
# the implementation of Geometry.get_cell_coordinates(trian::SkeletonGrid)
function Gridap.Arrays.return_cache(k::Reindex,i::VectorBlock{<:Vector{<:Integer}})
  T=Gridap.Arrays.return_type(evaluate,k,i.array[1])
  Tc=typeof(Gridap.Arrays.return_cache(evaluate,k,i.array[1]))
  r=Vector{T}(undef,length(i.array))
  rc=Vector{Tc}(undef,length(i.array))
  rc[1]=Gridap.Arrays.return_cache(evaluate,k,i.array[1])
  for j=2:length(i.array)
    rc[j]=Gridap.Arrays.return_cache(evaluate,k,i.array[j])
  end
  (Gridap.Fields.ArrayBlock(r,i.touched),Gridap.Fields.ArrayBlock(rc,i.touched))
end
function Gridap.Arrays.evaluate!(cache,k::Reindex,i::VectorBlock{<:Vector{<:Integer}})
  (r,rc)=cache
  Gridap.Helpers.@check length(r.array) == length(i.array)
  for j=1:length(i.array)
    r.array[j]=Gridap.Arrays.evaluate!(rc.array[j],k,i.array[j])
  end
  r
end
function Gridap.Arrays.evaluate(k::Reindex,i::VectorBlock{<:Vector{<:Integer}})
  cache=return_cache(k,i)
  evaluate!(cache,k,i)
end
# Default return_type fails because one(::VectorBlock{<:Vector{<:Integer}}) is NOT defined
# This justifies why I had to define the following function
function Gridap.Arrays.return_type(k::Reindex,i::VectorBlock{<:Vector{<:Integer}})
  TR=eltype(k.values)
  VectorBlock{Vector{TR}}
end
# Default return_values fails because one(::VectorBlock{<:Vector{<:Integer}}) is NOT defined
# This justifies why I had to define the following function
function Gridap.Arrays.return_value(k::Reindex,i::VectorBlock{<:Vector{<:Integer}})
  evaluate(k,i)
end

function Geometry.get_cell_map(trian::SkeletonGrid)
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_shapefuns = get_cell_shapefuns(trian)
  lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)
end

struct SkeletonTriangulation{Dc,Dp,A,B,C} <: Triangulation{Dc,Dp}
  model::A
  grid::B
  sign_flip::C
  glue::Gridap.Geometry.FaceToCellGlue
  function SkeletonTriangulation(model::DiscreteModel)
    A     = typeof(model)
    D     = num_cell_dims(model)
    mgrid = get_grid(model)
    fgrid = Grid(ReferenceFE{D-1},model)
    glue =  Gridap.Geometry.FaceToCellGlue(get_grid_topology(model),
                                           mgrid,
                                           fgrid,
                                           collect(1:num_facets(model)),
                                           Fill(Int8(1),num_facets(model)))
    sgrid = SkeletonGrid(mgrid)
    B = typeof(sgrid)

    # Generate sign_flip
    # TO-DO: here I am reusing the machinery for global RT FE spaces.
    #        Sure there is a way to decouple this from global RT FE spaces.
    function _get_sign_flip(model)
      basis,reffe_args,reffe_kwargs = ReferenceFE(raviart_thomas,Float64,0)
      cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
      Gridap.FESpaces.get_sign_flip(model,cell_reffe)
    end
    sign_flip=_get_sign_flip(model)

    C = typeof(sign_flip)
    new{D-1,D,A,B,C}(model,sgrid,sign_flip,glue)
  end
end

Skeleton(args...) = SkeletonTriangulation(args...)

Geometry.get_background_model(a::SkeletonTriangulation) = a.model
Geometry.get_grid(a::SkeletonTriangulation) = a.grid

struct SkeletonGlue{A,B}
  trian::SkeletonTriangulation
  tcell_lface_mface::A
  tcell_lface_mface_map::B
end

function Geometry.is_change_possible(sglue::FaceToFaceGlue,tglue::SkeletonGlue)
  true
end

function Geometry.is_change_possible(sglue::SkeletonGlue,tglue::SkeletonGlue)
  true
end

function Geometry.is_change_possible(
     strian::BodyFittedTriangulation{1},
     ttrian::BodyFittedTriangulation{2})
  get_background_model(strian) === get_background_model(ttrian)
end

function Geometry.is_change_possible(
    strian::BodyFittedTriangulation{2},
    ttrian::BodyFittedTriangulation{1})
  get_background_model(strian) === get_background_model(ttrian)
end

function Geometry.is_change_possible(
    strian::BodyFittedTriangulation{3},
    ttrian::BodyFittedTriangulation{2})
  get_background_model(strian) === get_background_model(ttrian)
end

function Geometry.is_change_possible(
    strian::BodyFittedTriangulation{2},
    ttrian::BodyFittedTriangulation{3})
  get_background_model(strian) === get_background_model(ttrian)
end

# In my view this is a little bit dirty. Transforming, e.g., the dof ids
# from a cell triangulation to a facet triangulation is not a well-defined
# operation. You loose information along the way, e.g., the interior DoFs.
# I had to define it because I defined is_change_possible(...) for this combination
# of triangulations. In turn, I needed to define is_change_possible(...)
# as I need to perform change_domain on fields defined at different triangulations,
# as per-required by Hybridizable methods
function Gridap.FESpaces.get_cell_fe_data(
  fun,
  sface_to_data,
  sglue::FaceToFaceGlue,
  tglue::Nothing)
  sface_to_data
end

function Geometry.best_target(a::BodyFittedTriangulation{Dca},
                              b::BodyFittedTriangulation{Dcb}) where {Dca,Dcb}
  @assert Dca==Dcb-1 || Dca-1==Dcb
  @assert get_background_model(a)===get_background_model(b)
  Skeleton(get_background_model(a))
end

# TO-DO: dirty. I cannot check whether a===b, as a and b might be created from scratch
#               along the process
function Geometry.best_target(a::SkeletonTriangulation{Dc},
                              b::SkeletonTriangulation{Dc}) where {Dc}
  a
end

function CellData.change_domain_ref_ref(a::CellField,
                                        ttrian::SkeletonTriangulation,
                                        sglue::SkeletonGlue,tglue::SkeletonGlue)
  a
end

function Geometry.get_glue(trian::SkeletonTriangulation{D},::Val{D}) where D
  model = get_background_model(trian)
  topo = get_grid_topology(model)
  cell_lface_face = get_faces(topo,D+1,D)
  pgrid = trian.grid.parent
  ctype_creffe = get_reffes(pgrid)
  ctype_lface_map = map(ctype_creffe) do reffe
    poly = get_polytope(reffe)
    fill(GenericField(identity),num_faces(poly,D))
  end
  cell_ctype = get_cell_type(pgrid)
  cell_lface_map = expand_cell_data(ctype_lface_map,cell_ctype)
  SkeletonGlue(trian,cell_lface_face,cell_lface_map)
end

function Geometry.get_glue(trian::SkeletonTriangulation{d},::Val{D}) where {d,D}
  if d+1 != D
    return nothing
  end
  pgrid = trian.grid.parent
  ctype_reffe = get_reffes(pgrid)
  cell_ctype = get_cell_type(pgrid)
  ncells = length(cell_ctype)
  cell_cell = IdentityVector(ncells)
  ctype_nlfaces = map(ctype_reffe) do reffe
    poly = get_polytope(reffe)
    num_faces(poly,d)
  end
  # Avoid allocations here
  tcell_lface_mface = lazy_map(cell_ctype,cell_cell) do ctype, cell
    nlfaces = ctype_nlfaces[ctype]
    fill(cell,nlfaces)
  end
  tcell_lface_mface_map = _setup_tcell_lface_mface_map(d,trian.model,trian.glue)
  SkeletonGlue(trian,tcell_lface_mface,tcell_lface_mface_map)
end

function CellData.change_domain_ref_ref(
  a::CellField,ttrian::SkeletonTriangulation,sglue::FaceToFaceGlue,tglue::SkeletonGlue)

  D = num_cell_dims(ttrian.model)
  strian = get_triangulation(a)

  @notimplementedif !(isa(strian,Triangulation{D-1,D}) || isa(strian,Triangulation{D,D}))

  if isa(strian,Triangulation{D,D})
    b=_restrict_to_cell_boundary_cell_fe_basis(ttrian.model,
                                                ttrian.glue,
                                                tglue.tcell_lface_mface_map,
                                                a)
  elseif isa(strian,Triangulation{D-1,D})
    b=_restrict_to_cell_boundary_facet_fe_basis(ttrian.model,ttrian.glue,a)
  end
  CellData.similar_cell_field(a,b,ttrian,ReferenceDomain())
end

function _restrict_to_cell_boundary_cell_fe_basis(model,
                                                  glue,
                                                  tface_to_mface_map,
                                                  cell_fe_basis::Gridap.CellData.CellField)
  D = num_cell_dims(model)
  Gridap.Helpers.@check isa(get_triangulation(cell_fe_basis),Triangulation{D,D})
  cell_a_q = _transform_cell_to_cell_lface_array(glue,
         Gridap.CellData.get_data(cell_fe_basis);
         add_naive_innermost_block_level=true)
  lazy_map(Broadcasting(∘),cell_a_q,tface_to_mface_map)
end

function _transform_cell_to_cell_lface_array(glue,
                                            cell_array::Fill;
                                            add_naive_innermost_block_level=false)
    d = Gridap.Arrays.CompressedArray([cell_array.value,],Fill(1,length(cell_array)))
    _transform_cell_to_cell_lface_array(glue,d;add_naive_innermost_block_level=add_naive_innermost_block_level)
end


function _transform_cell_to_cell_lface_array(glue,
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

function _transform_cell_to_cell_lface_array(
  glue,
  cell_array::Gridap.Arrays.LazyArray{<:Fill{typeof(transpose)}};
  add_naive_innermost_block_level=false)
  Gridap.Helpers.@check typeof(cell_array.args[1]) <: Fill
  cell_array_fill = Fill(evaluate(transpose,cell_array.args[1].value),length(cell_array))
  _transform_cell_to_cell_lface_array(glue,cell_array_fill; add_naive_innermost_block_level=add_naive_innermost_block_level)
end

function _restrict_to_cell_boundary_facet_fe_basis(model,
                                                   glue,
                                                   facet_fe_basis::Gridap.CellData.CellField)

  D = num_cell_dims(model)
  Gridap.Helpers.@check isa(get_triangulation(facet_fe_basis),Triangulation{D-1,D})

  _transform_face_to_cell_lface_expanded_array(
    glue,
    Gridap.CellData.get_data(facet_fe_basis))
end

function _transform_face_to_cell_lface_expanded_array(
  glue,
  face_array::Gridap.Arrays.LazyArray{<:Fill{typeof(transpose)}})
  Gridap.Helpers.@check typeof(face_array.args[1]) <: Gridap.Arrays.CompressedArray

  T = Gridap.Arrays.return_type(transpose,face_array.args[1].values[1])
  v = Vector{T}(undef,length(face_array.args[1].values))
  for i=1:length(face_array.args[1].values)
    v[i]=evaluate(transpose,face_array.args[1].values[i])
  end
  face_array_compressed=Gridap.Arrays.CompressedArray(v,face_array.args[1].ptrs)
  _transform_face_to_cell_lface_expanded_array(glue,face_array_compressed)
end


function _transform_face_to_cell_lface_expanded_array(glue,
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

function _transform_face_to_cell_lface_expanded_array(glue,
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

function _get_block_layout(fields_array::AbstractArray{<:AbstractArray{<:Gridap.Fields.Field}})
  Fill((1,1),length(fields_array))
end

function _get_block_layout(fields_array::AbstractArray{<:Gridap.Fields.ArrayBlock})
  lazy_map(x->((size(x),findall(x.touched))),fields_array)
end


function CellData.change_domain_phys_phys(
  a::CellField,ttrian::SkeletonTriangulation,sglue::FaceToFaceGlue,tglue::SkeletonGlue)
  sface_to_field = get_data(a)
  mface_to_sface = sglue.mface_to_tface
  tcell_lface_mface = tglue.tcell_lface_mface
  mface_to_field = extend(sface_to_field,mface_to_sface)
  # TODO this can be optimized
  tface_to_field = lazy_map(tcell_lface_mface) do lface_mface
    mface_to_field[lface_mface]
  end
  CellData.similar_cell_field(a,tface_to_field,ttrian,PhysicalDomain())
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell.
"""
function get_cell_normal_vector(s::SkeletonTriangulation)
  cell_lface_normal=_get_cell_normal_vector(s.model, s.glue, _cell_lface_to_nref)
  GenericCellField(cell_lface_normal,s,ReferenceDomain())
end

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
  SkeletonCompressedVector(ctype_lface_pindex_to_nref,glue)
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell owner of the facet.
"""
function get_cell_owner_normal_vector(s::SkeletonTriangulation)
  cell_owner_lface_normal=_get_cell_normal_vector(
     s.model,
     s.glue,
     _cell_lface_to_owner_nref,
     s.sign_flip)
  GenericCellField(cell_owner_lface_normal,s,ReferenceDomain())
end

function _cell_lface_to_owner_nref(args...)
  model,glue,sign_flip = args
  cell_lface_to_nref=_cell_lface_to_nref(model,glue)
  SkeletonOwnerNref(cell_lface_to_nref,sign_flip)
end

function _get_cell_normal_vector(model,glue,cell_lface_to_nref::Function,sign_flip=nothing)
  cell_grid = get_grid(model)

  cell_lface_to_nref = cell_lface_to_nref(model,glue,sign_flip)
  cell_lface_s_nref = lazy_map(Gridap.Fields.constant_field,cell_lface_to_nref)

  # Inverse of the Jacobian transpose
  cell_q_x = get_cell_map(cell_grid)
  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(Gridap.Fields.pinvJt),cell_q_Jt)
  cell_lface_q_invJt = _transform_cell_to_cell_lface_array(glue, cell_q_invJt)

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
  cell_lface_to_q_vertex_coords = SkeletonCompressedVector(
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

function Gridap.ReferenceFEs.expand_cell_data(
  type_to_data,
  cell_to_type::CompressedArray{T,1}) where T <: AbstractVector{<:Integer}

  TD=eltype(type_to_data)
  ctype_to_vector_block=
    Vector{Gridap.Fields.VectorBlock{TD}}(undef,length(cell_to_type.data))
  for ctype=1:length(cell_to_type.data)
     num_facets=length(cell_to_type.data[ctype])
     v=Vector{TD}(undef,num_facets)
     t=Vector{Bool}(undef,num_facets)
     t.=true
     for lface=1:num_facets
      ftype=cell_to_type.data[ctype][lface]
      v[lface]=type_to_data[ftype]
     end
     ctype_to_vector_block[ctype]=Gridap.Fields.ArrayBlock(v,t)
  end
  Gridap.Arrays.CompressedArray(ctype_to_vector_block,cell_to_ctype.ptrs)
end


function Gridap.ReferenceFEs.expand_cell_data(
  type_to_data,
  cell_to_type::Fill{T,1}) where T <: AbstractVector{<:Integer}

  TD=eltype(type_to_data)
  num_facets=length(cell_to_type.value)
  v=Vector{TD}(undef,num_facets)
  t=Vector{Bool}(undef,num_facets)
  t.=true
  for lface=1:num_facets
    ftype=cell_to_type.value[lface]
    v[lface]=type_to_data[ftype]
  end
  Fill(Gridap.Fields.ArrayBlock(v,t),length(cell_to_type))
end


# TODO this needs to be optimized (very important)
function Arrays.evaluate!(
  cache,
  f::AbstractVector,
  x::AbstractVector{<:AbstractVector{<:Point}})
  @check length(f) == length(x)
  evaluate.(f,x)
end

function Visualization.visualization_data(
  a::SkeletonTriangulation,
  filebase::AbstractString;
  offset=0,
  cellfields=Dict())

  model = get_background_model(a)
  grid = get_grid(a)
  D = num_cell_dims(model)
  node_coord = get_node_coordinates(grid)
  cell_lface_lfnode_node = get_cell_node_ids(grid)
  ncells = length(cell_lface_lfnode_node)
  parent = grid.parent
  cell_nodes = get_cell_node_ids(parent)

  P = eltype(node_coord)
  fnode_coord = P[]
  face_fnodes = Vector{Int32}[]
  fnode = Int32(0)
  for cell in 1:ncells
    lnode_node = cell_nodes[cell] # Allocation here
    lnode_coord = node_coord[lnode_node]
    Xm = sum(lnode_coord) / length(lnode_coord)
    lface_lfnode_node = cell_lface_lfnode_node[cell] # Allocation here
    for lfnode_node in lface_lfnode_node
      fnodes = Int32[] # Allocation here
      for node in lfnode_node
         Xf = node_coord[node]
         coord = Xf + offset*(Xm-Xf)
         push!(fnode_coord,coord) # Allocation here
         fnode += Int32(1)
         push!(fnodes,fnode) # Allocation here
      end
      push!(face_fnodes,fnodes) # Allocation here
    end
  end

  function compute_pdata(f)
    x = get_cell_points(a)
    cell_lface_node_val = f(x)
    T = eltype(eltype(eltype(cell_lface_node_val)))
    vals = zeros(T,length(fnode_coord))
    i = 0
    # To be optimized
    for lface_node_val in cell_lface_node_val
      for node_val in lface_node_val
        for val in node_val
          i += 1
          vals[i] = val
        end
      end
    end
    vals
  end

  pdata = Dict()
  for (k,v) in cellfields
    pdata[k] = compute_pdata(v)
  end

  ftype_freffe = get_reffes(grid)
  @notimplementedif length(ftype_freffe) != 1
  freffe = first(ftype_freffe)
  freffes = [freffe,]
  nfaces = length(face_fnodes)
  face_ftype = ones(Int8,nfaces)

  fgrid = UnstructuredGrid(
    fnode_coord,Table(face_fnodes),freffes,face_ftype)
  (VisualizationData(fgrid,filebase;nodaldata=pdata),)
end


function Gridap.CellData.change_domain(
  a::Gridap.MultiField.MultiFieldFEBasisComponent,
  ttrian::SkeletonTriangulation,
  tdomain::DomainStyle)
  cf_cell_basis=GenericCellField(a.cell_basis,get_triangulation(a),DomainStyle(a))
  b=change_domain(cf_cell_basis,ttrian,tdomain)
  Gridap.MultiField.MultiFieldFEBasisComponent(Gridap.CellData.get_data(b),b,a.fieldid,a.nfields)
end
