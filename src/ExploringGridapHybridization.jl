module ExploringGridapHybridization
using Base: _Set, ArithmeticRounds
using Gridap
using FillArrays
using LinearAlgebra
include("GridapOverloads.jl")
include("CellBoundary.jl")
include("DensifyInnerMostBlockLevel.jl")
include("StaticCondensationMap.jl")
include("BackwardStaticCondensationMap.jl")


export DensifyInnerMostBlockLevelMap
export StaticCondensationMap
export BackwardStaticCondensationMap

export CellBoundary
export get_cell_owner_normal_vector
export get_cell_normal_vector
export quadrature_evaluation_points_and_weights
export restrict_to_cell_boundary
export restrict_facet_dof_ids_to_cell_boundary

export integrate_vh_mult_uh_low_level
export integrate_vh_cdot_n_mult_lh_low_level
export integrate_mh_mult_uh_cdot_n_low_level

end # module
