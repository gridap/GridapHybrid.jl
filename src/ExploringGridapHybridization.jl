module ExploringGridapHybridization
using Base: _Set, ArithmeticRounds
using Gridap
using FillArrays
using LinearAlgebra

include("GridapAPIExtensions.jl")
include("GridapTmpModifications.jl")

using Gridap.Fields
include("StaticCondensationMap.jl")
include("BackwardStaticCondensationMap.jl")
include("RestrictArrayBlockMap.jl")
include("AddNaiveInnerMostBlockLevelMap.jl")
include("Scalar2ArrayBlockMap.jl")
include("SumFacetsMap.jl")

export StaticCondensationMap
export BackwardStaticCondensationMap
export RestrictArrayBlockMap
export AddNaiveInnerMostBlockLevelMap
export Scalar2ArrayBlockMap
export SumFacetsMap

export get_cell_owner_normal_vector
export get_cell_normal_vector

include("SkeletonArrays.jl")
include("Skeleton.jl")
export Skeleton

include("HybridAffineFEOperators.jl")
export HybridAffineFEOperator

end # module
