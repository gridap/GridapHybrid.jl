using Gridap
using FillArrays
using Gridap.Geometry

struct CellBoundary{M<:DiscreteModel,BTP<:Triangulation,BTM<:Triangulation}
  model::M
  boundary_trian_plus::BTP
  boundary_trian_minus::BTM
  sign_flip

  function CellBoundary(model::M) where M<:DiscreteModel

     function _get_bgface_to_lcell_btm(model)
       bgface_to_lcell_btm  = Vector{Int8}(undef,num_facets(model))
       bgface_to_lcell_btm     .= 2

       Dc=num_cell_dims(model)
       Df=Dc-1
       topo = Gridap.Geometry.get_grid_topology(model)
       bgface_to_cells = Gridap.Geometry.get_faces(topo,Df,Dc)

       c=array_cache(bgface_to_cells)
       for i=1:length(bgface_to_cells)
         cells_around=getindex!(c,bgface_to_cells,i)
         if length(cells_around)==1
          bgface_to_lcell_btm[i] = 1
         end
       end
       bgface_to_lcell_btm
     end

     # TO-DO: here I am reusing the machinery for global RT FE spaces.
     #        Sure there is a way to decouple this from global RT FE spaces.
     function _get_sign_flip(model)
       basis,reffe_args,reffe_kwargs = ReferenceFE(raviart_thomas,Float64,0)
       cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
       Gridap.FESpaces.get_sign_flip(model,cell_reffe)
     end


    face_to_bgface       = collect(1:num_facets(model))
    btp=BoundaryTriangulation(model,
                              face_to_bgface,
                              Fill(Int8(1),num_facets(model)))
    btm=BoundaryTriangulation(model,
                              face_to_bgface,
                              _get_bgface_to_lcell_btm(model))



    TBTP=typeof(btp)
    TBTM=typeof(btm)
    new{M,TBTP,TBTM}(model,btp,btm,_get_sign_flip(model))
  end
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell.
"""
function get_cell_normal_vector(cb::CellBoundary)
    np=get_normal_vector(cb.boundary_trian_plus)
    nm=get_normal_vector(cb.boundary_trian_minus)
    signed_cell_wise_facets_ids = _get_signed_cell_wise_facets(cb)
    np_array=Gridap.CellData.get_data(np)
    nm_array=Gridap.CellData.get_data(nm)
    num_cell_facets = _get_num_facets(cb)
    per_local_facet_normal_cell_wise_arrays =
       _get_per_local_facet_cell_wise_arrays(np_array,
                                             nm_array,
                                             signed_cell_wise_facets_ids,
                                             num_cell_facets)
    _block_arrays(per_local_facet_normal_cell_wise_arrays...)
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell owner of the facet.
"""
function get_cell_owner_normal_vector(cb::CellBoundary)
  np=get_normal_vector(cb.boundary_trian_plus)
  np_array=Gridap.CellData.get_data(np)
  num_cell_facets = _get_num_facets(cb)
  cell_wise_facets = _get_cell_wise_facets(cb)
  per_local_facet_normal_cell_wise_arrays =
      _get_per_local_facet_cell_wise_arrays(np_array,
                                            cell_wise_facets,
                                            num_cell_facets)
  _block_arrays(per_local_facet_normal_cell_wise_arrays...)
end

function quadrature_evaluation_points_and_weights(cell_boundary::CellBoundary, degree::Integer)
  model = cell_boundary.model
  D = num_cell_dims(model)
  p = lazy_map(Gridap.ReferenceFEs.get_polytope,get_reffes(model))
  Gridap.Helpers.@check length(p) == 1
  p  = p[1]
  pf = Gridap.ReferenceFEs.Polytope{D-1}(p,1)
  qf = Gridap.ReferenceFEs.Quadrature(pf,degree)

  xf = Gridap.ReferenceFEs.get_coordinates(qf)
  wf = Gridap.ReferenceFEs.get_weights(qf)

  num_cell_facets = _get_num_facets(cell_boundary)
  xf_array_block=_block_arrays(Fill(Fill(xf,num_cells(model)),num_cell_facets)...)
  wf_array_block=_block_arrays(Fill(Fill(wf,num_cells(model)),num_cell_facets)...)

  (xf_array_block,wf_array_block)
end

# This function plays a similar role of the change_domain function in Gridap.
# Except it does not, by now, deal with Triangulation and DomainStyle instances.
# (TO-THINK)
function restrict_to_cell_boundary(cb::CellBoundary,
                                   fe_basis::Gridap.FESpaces.FEBasis)
  model = cb.model
  D = num_cell_dims(model)
  trian = get_triangulation(fe_basis)
  if isa(trian,Triangulation{D-1,D})
    _restrict_to_cell_boundary_facet_fe_basis(cb,fe_basis)
  elseif isa(trian,Triangulation{D,D})
    _restrict_to_cell_boundary_cell_fe_basis(cb,fe_basis)
  end
end

function restrict_to_cell_boundary(cb::CellBoundary, cell_field::Gridap.CellData.CellField)
  model = cb.model
  D = num_cell_dims(model)
  trian = get_triangulation(cell_field)
  if isa(trian,Triangulation{D-1,D})
    # _restrict_to_cell_boundary_facet_fe_basis(cb,cell_field)
    @assert false
  elseif isa(trian,Triangulation{D,D})
    _restrict_to_cell_boundary_cell_fe_basis(cb,cell_field)
  end
end


function _restrict_to_cell_boundary_facet_fe_basis(cb::CellBoundary,
                                                   facet_fe_basis::Gridap.FESpaces.FEBasis;
                                                   facet_last_block_level=true)

  # TO-THINK:
  #     1) How to deal with CellField objects which are NOT FEBasis objects?
  #     2) How to deal with differential operators applied to FEBasis/CellField objects?
  D = num_cell_dims(cb.model)
  Gridap.Helpers.@check isa(get_triangulation(facet_fe_basis),Triangulation{D-1,D})

  num_cell_facets  = _get_num_facets(cb)
  cell_wise_facets = _get_cell_wise_facets(cb)

  fta=Gridap.CellData.get_data(facet_fe_basis)
  fe_basis_restricted_to_lfacet_cell_wise =
    _get_per_local_facet_cell_wise_arrays(fta,
                                          cell_wise_facets,
                                          num_cell_facets)
  fe_basis_restricted_to_lfacet_cell_wise_expanded =
         [ _expand_facet_lid_fields_to_num_facets_blocks(
                            Gridap.FESpaces.BasisStyle(facet_fe_basis),
                            fe_basis_restricted_to_lfacet_cell_wise[facet_lid],
                            num_cell_facets,
                            facet_lid) for facet_lid=1:num_cell_facets ]

  # Cell-wise array of VectorBlock entries with as many entries as facets.
  # For each cell, and local facet, we get the cell FE Basis restricted to that facet.
  _block_arrays(fe_basis_restricted_to_lfacet_cell_wise_expanded...)
end

function _get_block_layout(fields_array::AbstractArray{<:AbstractArray{<:Gridap.Fields.Field}})
  Fill((1,1),length(fields_array))
end

function _get_block_layout(fields_array::AbstractArray{<:Gridap.Fields.ArrayBlock})
  lazy_map(x->((size(x),findall(x.touched))),fields_array)
end

function _restrict_to_cell_boundary_cell_fe_basis(cb::CellBoundary,
                                                  cell_fe_basis::Union{Gridap.FESpaces.FEBasis,Gridap.CellData.CellField})
  # TO-THINK:
  #     1) How to deal with CellField objects which are NOT FEBasis objects?
  #     2) How to deal with differential operators applied to FEBasis objects?
  D = num_cell_dims(cb.model)
  Gridap.Helpers.@check isa(get_triangulation(cell_fe_basis),Triangulation{D,D})

  num_cell_facets = _get_num_facets(cb)

  cfp=Gridap.CellData.change_domain(cell_fe_basis,
                                    cb.boundary_trian_plus,
                                    DomainStyle(cell_fe_basis))
  cfm=Gridap.CellData.change_domain(cell_fe_basis,
                                    cb.boundary_trian_minus,
                                    DomainStyle(cell_fe_basis))
  cfp_array=Gridap.CellData.get_data(cfp)
  cfm_array=Gridap.CellData.get_data(cfm)
  signed_cell_wise_facets_ids = _get_signed_cell_wise_facets(cb)

  # Array with as many cell-wise arrays as facets. For each facet, we get a
  # cell-wise array with the cell FE Basis restricted to that facet.
  per_local_facet_fe_basis_cell_wise_arrays =
     _get_per_local_facet_cell_wise_arrays(cfp_array,
                                           cfm_array,
                                           signed_cell_wise_facets_ids,
                                           num_cell_facets)

  # Cell-wise array of VectorBlock entries with as many entries as facets.
  # For each cell, and local facet, we get the cell FE Basis restricted to that facet.
  _block_arrays(per_local_facet_fe_basis_cell_wise_arrays...)
end

function _expand_facet_lid_fields_to_num_facets_blocks(::Gridap.FESpaces.TrialBasis,
                                                       fe_basis_restricted_to_lfacet_cell_wise,
                                                       num_cell_facets,
                                                       facet_lid)
    bl=_get_block_layout(fe_basis_restricted_to_lfacet_cell_wise)[1]
    lazy_map(Gridap.Fields.BlockMap(bl[1],bl[2]),
             lazy_map(Gridap.Fields.BlockMap((1,num_cell_facets),[CartesianIndex((1,facet_lid))]),
                      lazy_map(x->x[findfirst(x.touched)], fe_basis_restricted_to_lfacet_cell_wise)))
end

function _expand_facet_lid_fields_to_num_facets_blocks(::Gridap.FESpaces.TestBasis,
                                                       fe_basis_restricted_to_lfacet_cell_wise,
                                                       num_cell_facets,
                                                       facet_lid)
    bl=_get_block_layout(fe_basis_restricted_to_lfacet_cell_wise)[1]
    lazy_map(Gridap.Fields.BlockMap(bl[1][1],bl[2]),
             lazy_map(Gridap.Fields.BlockMap(num_cell_facets,facet_lid),
                      lazy_map(x->x[findfirst(x.touched)], fe_basis_restricted_to_lfacet_cell_wise)))
end

function _get_num_facets(cb::CellBoundary)
  model = cb.model
  p = lazy_map(Gridap.ReferenceFEs.get_polytope,get_reffes(model))
  Gridap.Helpers.@check length(p) == 1
  num_facets(p[1])
end

function _get_cell_wise_facets(cb::CellBoundary)
  model = cb.model
  gtopo = get_grid_topology(model)
  D     = num_cell_dims(model)
  Gridap.Geometry.get_faces(gtopo, D, D-1)
end

function _get_cells_around_facets(cb::CellBoundary)
  model = cb.model
  gtopo = get_grid_topology(model)
  D     = num_cell_dims(model)
  Gridap.Geometry.get_faces(gtopo, D-1, D)
end

function _get_signed_cell_wise_facets(cb::CellBoundary)
  lazy_map(Broadcasting((y,x)-> y ? -x : x),
                        cb.sign_flip,
                        _get_cell_wise_facets(cb))
end

function Gridap.Geometry.get_cell_map(cb::CellBoundary)
  num_cell_facets = _get_num_facets(cb)
  cell_map_tp=get_cell_map(cb.boundary_trian_plus)
  cell_map_tm=get_cell_map(cb.boundary_trian_minus)
  signed_cell_wise_facets_ids = _get_signed_cell_wise_facets(cb)
  # Array with as many cell-wise arrays as facets. For each facet, we get a
  # cell-wise array with the cell FE Basis restricted to that facet.
  per_local_facet_maps_cell_wise_arrays =
      _get_per_local_facet_cell_wise_arrays(cell_map_tp,
                                            cell_map_tm,
                                            signed_cell_wise_facets_ids,
                                            num_cell_facets)
  _block_arrays(per_local_facet_maps_cell_wise_arrays...)
end

function _get_per_local_facet_cell_wise_arrays(tpa, # Facet triangulation + array
                                               tma, # Facet triangulation - array
                                               signed_cell_wise_facets_ids,
                                               num_cell_facets)
  [ _get_local_facet_cell_wise_array(tpa,
                                     tma,
                                     signed_cell_wise_facets_ids,
                                     facet_lid) for facet_lid=1:num_cell_facets ]
end

function _get_local_facet_cell_wise_array(tpa,
                                          tma,
                                          signed_cell_wise_facets_ids,
                                          facet_lid)
  Gridap.Helpers.@check length(tpa) == length(tma)
  tpa_tma=Gridap.Arrays.AppendedArray(tpa,tma)
  function transform_signed_index_to_appended_index(x)
      signed_index=x[facet_lid]
      Gridap.Helpers.@check signed_index != 0
      if (signed_index>0)
        return signed_index
      else
        return length(tpa)+abs(signed_index)
      end
  end
  indices=lazy_map(transform_signed_index_to_appended_index,
                   signed_cell_wise_facets_ids)
  lazy_map(Gridap.Arrays.Reindex(tpa_tma),indices)
end

function _get_per_local_facet_cell_wise_arrays(fta, # Facet triangulation array
                                               cell_wise_facets_ids,
                                               num_cell_facets)
  [ _get_local_facet_cell_wise_array(fta,
                                     cell_wise_facets_ids,
                                     facet_lid) for facet_lid=1:num_cell_facets ]
end

function _get_local_facet_cell_wise_array(fta,
                                          cell_wise_facets_ids,
                                          facet_lid)
  indices=lazy_map((x)->x[facet_lid],cell_wise_facets_ids)
  lazy_map(Gridap.Arrays.Reindex(fta),indices)
end

function _block_arrays(a::AbstractArray...)
  num_blocks = length(a)
  lazy_map(Gridap.Fields.BlockMap(num_blocks,collect(1:num_blocks)),a...)
end

# This function will eventually play the role of its counterpart in Gridap
# function integrate(f::CellField,quad::CellQuadrature)
# TO-THINK: mh,uh currently as LazyArray{...}
#    Not sure if we are loosing generality by constraining a to be of type
#    LazyArray{...}. I need to restrict the cell-wise block array to each
#    individual block, and with a LazyArray{...} this is very efficient as
#    the array is already restricted to each block in the a.args member variable.

#∫( mh*(uh⋅n) )*dK
function integrate_mh_mult_uh_cdot_n_low_level(cb::CellBoundary,
  mh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  uh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
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

  sum_facets=_sum_facets(cb,mhx_mult_uhx_cdot_nx,jx,w)
  lazy_map(DensifyInnerMostBlockLevelMap(),sum_facets)
end

function _sum_facets(cb,vx,jx,w)
 result=lazy_map(Broadcasting(+),
                 _set_up_integrate_block(vx,w,jx,1),
                 _set_up_integrate_block(vx,w,jx,2))
 for i=3:_get_num_facets(cb)
   result = lazy_map(Broadcasting(+),result,
                 _set_up_integrate_block(vx,w,jx,i))
 end
 result
end

#∫( (vh⋅n)*lh )*dK
function integrate_vh_cdot_n_mult_lh_low_level(
  cb::CellBoundary,
  vh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  lh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
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

  lazy_map(DensifyInnerMostBlockLevelMap(),_sum_facets(cb,vhx_cdot_nx_mult_lhx,jx,w))
end

# ∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
function integrate_qh_mult_uh_cdot_n_plus_stab_low_level(
  cb::CellBoundary,
  qh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  uh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  τ ::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  ph::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  lh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
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
                _sum_facets(cb,qhx_mult_uhx_cdot_nx,jx,w),
                _sum_facets(cb,qhx_mult_τx_mult_phx_mult_n_cdot_no,jx,w))

  arg2=lazy_map(DensifyInnerMostBlockLevelMap(),
                _sum_facets(cb,qhx_mult_τx_mult_lhx_mult_n_cdot_no,jx,w))

  lazy_map(Broadcasting(-),arg1,arg2)
end

# ∫(mh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
function integrate_mh_mult_uh_cdot_n_plus_stab_low_level(
  cb::CellBoundary,
  mh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  uh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  τ::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  ph::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
  lh::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
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

  arg1=lazy_map(DensifyInnerMostBlockLevelMap(),
               _sum_facets(cb,mhx_mult_uhx_cdot_nx,jx,w))
  arg2=lazy_map(DensifyInnerMostBlockLevelMap(),
               _sum_facets(cb,mhx_mult_τx_mult_phx_mult_n_cdot_no,jx,w))
  arg3=lazy_map(DensifyInnerMostBlockLevelMap(),
               _sum_facets(cb,mhx_mult_τx_mult_lhx_mult_n_cdot_no,jx,w))

  result=lazy_map(Broadcasting(+),arg1,arg2)
  result=lazy_map(Broadcasting(-),result,arg3)
end

function _set_up_integrate_block(intq,w,jq,block)
  lazy_map(Gridap.Fields.IntegrationMap(),
           _restrict_cell_array_block_to_block(intq,block),
           _restrict_cell_array_block_to_block(w,block),
           _restrict_cell_array_block_to_block(jq,block))
end

function restrict_facet_dof_ids_to_cell_boundary(cb::CellBoundary,facet_dof_ids)
  cell_wise_facets = _get_cell_wise_facets(cb)
  num_cell_facets = _get_num_facets(cb)

  facet_dof_ids_restricted_to_lfacet_cell_wise =
  _get_per_local_facet_cell_wise_arrays(facet_dof_ids,
                                        cell_wise_facets,
                                        num_cell_facets)

  # TO-DO: define LazyArray cache for mycat
  function mycat(x::Union{<:Vector{<:Integer},Integer}...)
    cat(x...,dims=1)
  end

  lazy_map(mycat,facet_dof_ids_restricted_to_lfacet_cell_wise...)
end


function _generate_glue_among_facet_and_cell_wise_dofs_arrays(cb::CellBoundary,facet_dof_ids)
  cells_around_facets=_get_cells_around_facets(cb)
  c1=array_cache(cells_around_facets)
  cell_wise_facets=_get_cell_wise_facets(cb)
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
       cb::CellBoundary,
       cell_dofs_array::AbstractVector{<:AbstractVector},
       facet_dof_ids)
  glue = _generate_glue_among_facet_and_cell_wise_dofs_arrays(cb,facet_dof_ids)
  k=ExtractFacetDofsFromCellDofs(cell_dofs_array)
  lazy_map(k,glue)
end
