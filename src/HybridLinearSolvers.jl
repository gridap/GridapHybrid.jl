struct HybridBackslashSymbolicSetup <: Gridap.Algebra.SymbolicSetup end

mutable struct HybridBackslashNumericalSetup{T<:HybridFEMatrix} <: Gridap.Algebra.NumericalSetup
  A::T
end

Gridap.Algebra.symbolic_setup(::BackslashSolver,mat::HybridFEMatrix) = HybridBackslashSymbolicSetup()

Gridap.Algebra.numerical_setup(::HybridBackslashSymbolicSetup,mat::HybridFEMatrix) = HybridBackslashNumericalSetup(mat)

function Gridap.Algebra.numerical_setup!(ns::HybridBackslashNumericalSetup, mat::HybridFEMatrix)
  ns.A = mat
end

function Gridap.Algebra.solve!(
  x::AbstractVector,ns::HybridBackslashNumericalSetup,b::HybridFEVector)

  Gridap.Helpers.@check ns.A.trial_hybridizable === b.trial_hybridizable
  Gridap.Helpers.@check ns.A.test_hybridizable  === b.test_hybridizable
  Gridap.Helpers.@check ns.A.trial_skeleton     === b.trial_skeleton
  Gridap.Helpers.@check ns.A.test_skeleton      === b.test_skeleton

  obiform, oliform = _merge_bulk_and_skeleton_contributions(ns.A.mat_contribs,b.vec_contribs)

  # Pair LHS and RHS terms associated to SkeletonTriangulation
  matvec,mat,vec=Gridap.FESpaces._pair_contribution_when_possible(obiform,oliform)

  # Add StaticCondensationMap to matvec terms
  matvec_statically_condensed=_add_static_condensation(matvec,ns.A.bulk_fields,ns.A.skeleton_fields)

  if (length(ns.A.skeleton_fields)!=1)
    matvec_block=_block_skeleton_system_contributions(matvec_statically_condensed,ns.A.test_skeleton)
  else
    matvec_block=matvec_statically_condensed
  end

  data = Gridap.FESpaces._collect_cell_matrix_and_vector(ns.A.trial_skeleton,
                                                         ns.A.test_skeleton,
                                                         matvec_block,
                                                         mat,
                                                         vec)

  A,b = assemble_matrix_and_vector(
    SparseMatrixAssembler(ns.A.trial_skeleton,ns.A.test_skeleton),data)
  x_skel = A\b

  free_dof_values=_compute_hybridizable_from_skeleton_free_dof_values(
                    x_skel,
                    ns.A.trial_hybridizable,
                    ns.A.test_hybridizable,
                    ns.A.trial_skeleton,
                    matvec,
                    ns.A.bulk_fields,
                    ns.A.skeleton_fields)
end
