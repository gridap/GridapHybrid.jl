module GridapHybrid
using Base: _Set, ArithmeticRounds
using Gridap
using FillArrays
using LinearAlgebra

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
include("HybridFEOperators.jl")
export HybridFEOperator

include("HybridLinearSolvers.jl")

include("GridapAPIExtensions.jl")
include("GridapTmpModifications.jl")

# HHO-methods specific components
include("OrthogonalBasisRefFEs.jl")
include("OrthogonalBasisFESpaces.jl")
export OrthogonalBasisRefFE
export orthogonal_basis

include("MonomialBasisRefFEs.jl")
export MonomialBasisRefFE
export monomial_basis

include("LocalFEOperators.jl")
export LocalFEOperator
export SingleValued
export MultiValued

end # module
