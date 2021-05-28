module ExploringGridapHybridization
using Gridap
using FillArrays
using LinearAlgebra
include("GridapOverloads.jl")
include("CellBoundary.jl")
include("DensifyInnerMostBlockLevel.jl")
include("StaticCondensationMap.jl")

export DensifyInnerMostBlockLevelMap
export StaticCondensationMap

export CellBoundary
export get_cell_owner_normal_vector
export get_cell_normal_vector
export quadrature_evaluation_points_and_weights
export restrict_to_cell_boundary

export integrate_vh_mult_uh_low_level
export integrate_vh_cdot_n_mult_lh_low_level
export integrate_mh_mult_uh_cdot_n_low_level

end # module
