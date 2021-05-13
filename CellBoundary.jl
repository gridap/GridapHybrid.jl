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


# TO-DO: discuss with @fverdugo
# I needed to overload this function in order to be able to build the LazyArray below
#    lazy_map(Broadcasting(Reindex(np_array)),cell_wise_facets_ids)

function Gridap.Arrays.return_cache(f::Broadcasting,x::Union{Number,AbstractArray{<:Number}}...)
  s = map(Gridap.Arrays._size,x)
  bs = Base.Broadcast.broadcast_shape(s...)
  T = Gridap.Arrays.return_type(f.f,map(Gridap.Arrays.testitem,x)...)
  N = length(bs)
  r = fill(Gridap.Arrays.return_value(f.f,map(Gridap.Arrays.testitem,x)...),bs)
  cache = Gridap.Arrays.CachedArray(r)
  Gridap.Arrays._prepare_cache!(cache,x...)
  cache
end


"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell.
"""
function get_cell_normal_vector(cell_boundary::CellBoundary)
    model = cell_boundary.model
    D = num_cell_dims(model)
    gtopo = get_grid_topology(model)
    # Extract composition among cells and facets
    cell_wise_facets_ids = Gridap.Geometry.get_faces(gtopo, D, D-1)
    np=get_normal_vector(cell_boundary.boundary_trian_plus)
    nm=get_normal_vector(cell_boundary.boundary_trian_minus)
    signed_cell_wise_facets_ids = lazy_map(Broadcasting((y,x)-> y ? -x : x),
                                           cell_boundary.sign_flip,
                                           cell_wise_facets_ids)
    np_array=Gridap.CellData.get_data(np)
    nm_array=Gridap.CellData.get_data(nm)
    lazy_map(Broadcasting(Gridap.Arrays.PosNegReindex(np_array,nm_array)),
             signed_cell_wise_facets_ids)
end

"""
Returns a cell-wise array which, for each cell, and each facet within the cell,
returns the unit outward normal to the boundary of the cell owner of the facet.
"""
function get_cell_owner_normal_vector(cell_boundary::CellBoundary)
  model = cell_boundary.model
  D = num_cell_dims(model)
  gtopo = get_grid_topology(model)
  # Extract composition among cells and facets
  cell_wise_facets_ids = Gridap.Geometry.get_faces(gtopo, D, D-1)
  np=get_normal_vector(cell_boundary.boundary_trian_plus)
  np_array=Gridap.CellData.get_data(np)
  lazy_map(Broadcasting(Reindex(np_array)),cell_wise_facets_ids)
end


domain = (0,1,0,1)
cells  = (1,2)
model  = CartesianDiscreteModel(domain,cells)
∂T     = CellBoundary(model)
nowner = get_cell_owner_normal_vector(∂T)
n      = get_cell_normal_vector(∂T)
