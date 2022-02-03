using Gridap.Helpers
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Arrays
using Gridap.Visualization
using FillArrays

export TempSkeleton

struct TempSkeletonGrid{Dc,Dp,P,A,B,C,D} <: Grid{Dc,Dp}
  parent::P
  node_coord::A
  cell_lface_nodes::B
  ftype_freffe::C
  cell_lface_ftype::D
  function TempSkeletonGrid(grid::Grid)
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

Geometry.get_node_coordinates(a::TempSkeletonGrid) = a.node_coord
Geometry.get_cell_node_ids(a::TempSkeletonGrid) = a.cell_lface_nodes
Geometry.get_reffes(a::TempSkeletonGrid) = a.ftype_freffe
Geometry.get_cell_type(a::TempSkeletonGrid) = a.cell_lface_ftype

# The following ones are a bit "hacky"

function Geometry.get_cell_coordinates(trian::TempSkeletonGrid)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_node_ids(trian)
  lazy_map(Reindex(node_to_coords),cell_to_nodes)
end

function Geometry.get_cell_map(trian::TempSkeletonGrid)
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_shapefuns = get_cell_shapefuns(trian)
  lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)
end

struct TempSkeletonTriangulation{Dc,Dp,A,B,C} <: Triangulation{Dc,Dp}
  model::A
  grid::B
  sign_flip::C
  glue::Gridap.Geometry.FaceToCellGlue
  function TempSkeletonTriangulation(model::DiscreteModel)
    A     = typeof(model)
    D     = num_cell_dims(model)
    mgrid = get_grid(model)
    fgrid = Grid(ReferenceFE{D-1},model)
    glue =  Gridap.Geometry.FaceToCellGlue(get_grid_topology(model),
                                           mgrid,
                                           fgrid,
                                           collect(1:num_facets(model)),
                                           Fill(Int8(1),num_facets(model)))
    sgrid = TempSkeletonGrid(mgrid)
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

TempSkeleton(args...) = TempSkeletonTriangulation(args...)

Geometry.get_background_model(a::TempSkeletonTriangulation) = a.model
Geometry.get_grid(a::TempSkeletonTriangulation) = a.grid

struct TempSkeletonGlue{A,B}
  trian::TempSkeletonTriangulation
  tcell_lface_mface::A
  tcell_lface_mface_map::B
end

function Geometry.is_change_possible(sglue::FaceToFaceGlue,tglue::TempSkeletonGlue)
  true
end

function Geometry.get_glue(trian::TempSkeletonTriangulation{D},::Val{D}) where D
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
  TempSkeletonGlue(trian,cell_lface_face,cell_lface_map)
end

function Geometry.get_glue(trian::TempSkeletonTriangulation{d},::Val{D}) where {d,D}
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
  TempSkeletonGlue(trian,tcell_lface_mface,tcell_lface_mface_map)
end

function CellData.change_domain_ref_ref(
  a::CellField,ttrian::TempSkeletonTriangulation,sglue::FaceToFaceGlue,tglue::TempSkeletonGlue)

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


function CellData.change_domain_phys_phys(
  a::CellField,ttrian::TempSkeletonTriangulation,sglue::FaceToFaceGlue,tglue::TempSkeletonGlue)
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

# TODO this needs to be optimized (very important)
function Arrays.evaluate!(
  cache,
  f::AbstractVector,
  x::AbstractVector{<:AbstractVector{<:Point}})
  @check length(f) == length(x)
  evaluate.(f,x)
end

function Visualization.visualization_data(
  a::TempSkeletonTriangulation,
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

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell.
"""
function get_cell_normal_vector(s::TempSkeletonTriangulation)
  cell_lface_normal=_get_cell_normal_vector(s.model, s.glue, _cell_lface_to_nref)
  GenericCellField(cell_lface_normal,s,ReferenceDomain())
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell owner of the facet.
"""
function get_cell_owner_normal_vector(s::TempSkeletonTriangulation)
  cell_owner_lface_normal=_get_cell_normal_vector(
     s.model,
     s.glue,
     _cell_lface_to_owner_nref,
     s.sign_flip)
  GenericCellField(cell_owner_lface_normal,s,ReferenceDomain())
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

function Gridap.CellData.CellQuadrature(trian::Triangulation,cell_quad,ids::DomainStyle)
  ctype_to_quad, cell_to_ctype = compress_cell_data(cell_quad)
  ctype_to_point = map(get_coordinates,ctype_to_quad)
  ctype_to_weight = map(get_weights,ctype_to_quad)
  cell_point = expand_cell_data(ctype_to_point,cell_to_ctype)
  cell_weight = expand_cell_data(ctype_to_weight,cell_to_ctype)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,ReferenceDomain(),ids)
end

function _get_f(v::VectorBlock,f::Function)
  n  = length(v.array)
  v1 = f(v.array[1])
  TD = typeof(v1)
  a  = Vector{TD}(undef,n)
  t  = v.touched
  for i=1:n
    a[i]=f(v.array[i])
  end
  Gridap.Fields.ArrayBlock(a,t)
end

function Gridap.Geometry.get_coordinates(v::VectorBlock)
  _get_f(v,Gridap.Geometry.get_coordinates)
end

function Gridap.Geometry.get_weights(v::VectorBlock)
  _get_f(v,Gridap.Geometry.get_weights)
end

function Gridap.CellData.change_domain(
  a::Gridap.MultiField.MultiFieldFEBasisComponent,
  ttrian::TempSkeletonTriangulation,
  tdomain::DomainStyle)
  cf_cell_basis=GenericCellField(a.cell_basis,get_triangulation(a),DomainStyle(a))
  b=change_domain(cf_cell_basis,ttrian,tdomain)
  Gridap.MultiField.MultiFieldFEBasisComponent(Gridap.CellData.get_data(b),b,a.fieldid,a.nfields)
end
