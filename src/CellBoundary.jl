using Gridap
using FillArrays
using Gridap.Geometry


function _get_block_layout(fields_array::AbstractArray{<:AbstractArray{<:Gridap.Fields.Field}})
  Fill((1,1),length(fields_array))
end

function _get_block_layout(fields_array::AbstractArray{<:Gridap.Fields.ArrayBlock})
  lazy_map(x->((size(x),findall(x.touched))),fields_array)
end


function _get_cell_wise_facets(cb)
  model = cb.model
  gtopo = get_grid_topology(model)
  D     = num_cell_dims(model)
  Gridap.Geometry.get_faces(gtopo, D, D-1)
end

function _get_cells_around_facets(cb)
  model = cb.model
  gtopo = get_grid_topology(model)
  D     = num_cell_dims(model)
  Gridap.Geometry.get_faces(gtopo, D-1, D)
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

  sum_facets=_sum_facets(cb,mhx_mult_uhx_cdot_nx,w,jx)
  lazy_map(DensifyInnerMostBlockLevelMap(),sum_facets)
end

function _sum_facets(cb::CellBoundary,vx,w,jx)
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

  lazy_map(DensifyInnerMostBlockLevelMap(),_sum_facets(cb,vhx_cdot_nx_mult_lhx,w,jx))
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
                _sum_facets(cb,qhx_mult_uhx_cdot_nx,w,jx),
                _sum_facets(cb,qhx_mult_τx_mult_phx_mult_n_cdot_no,w,jx))

  arg2=lazy_map(DensifyInnerMostBlockLevelMap(),
                _sum_facets(cb,qhx_mult_τx_mult_lhx_mult_n_cdot_no,w,jx))

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

  arg1=lazy_map(DensifyInnerMostBlockLevelMap(),
               _sum_facets(cb,mhx_mult_uhx_cdot_nx,w,jx))
  arg2=lazy_map(DensifyInnerMostBlockLevelMap(),
               _sum_facets(cb,mhx_mult_τx_mult_phx_mult_n_cdot_no,w,jx))
  arg3=lazy_map(DensifyInnerMostBlockLevelMap(),
               _sum_facets(cb,mhx_mult_τx_mult_lhx_mult_n_cdot_no,w,jx))

  result=lazy_map(Broadcasting(+),arg1,arg2)
  result=lazy_map(Broadcasting(-),result,arg3)
end

function _set_up_integrate_block(intq,w,jq,block)
  lazy_map(Gridap.Fields.IntegrationMap(),
           _restrict_cell_array_block_to_block(intq,block),
           _restrict_cell_array_block_to_block(w,block),
           _restrict_cell_array_block_to_block(jq,block))
end

function _generate_glue_among_facet_and_cell_wise_dofs_arrays(cb,facet_dof_ids)
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
       cb,
       cell_dofs_array::AbstractVector{<:AbstractVector},
       facet_dof_ids)
  glue = _generate_glue_among_facet_and_cell_wise_dofs_arrays(cb,facet_dof_ids)
  k=ExtractFacetDofsFromCellDofs(cell_dofs_array)
  lazy_map(k,glue)
end
