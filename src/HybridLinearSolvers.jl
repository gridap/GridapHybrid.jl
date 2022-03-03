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

  # Correction has to be zero at Dirichlet DoFs
  skeleton_fe_function = _generate_skeleton_fe_function(ns.A.trial_skeleton,x_skel)

  copy!(x,_compute_hybridizable_from_skeleton_free_dof_values(
                    skeleton_fe_function,
                    ns.A.trial_hybridizable,
                    ns.A.test_hybridizable,
                    ns.A.trial_skeleton,
                    matvec,
                    ns.A.bulk_fields,
                    ns.A.skeleton_fields))
  x
end

function _generate_skeleton_fe_function(fe::Gridap.FESpaces.SingleFieldFESpace, free_dofs)
  FEFunction(fe,free_dofs,zeros(num_dirichlet_dofs(fe)))
end

function _generate_skeleton_fe_function(fe::Gridap.MultiField.MultiFieldFESpace, free_dofs)
  blocks = Gridap.FESpaces.SingleFieldFEFunction[]
  for (field,U) in enumerate(fe.spaces)
    free_dofs_i = Gridap.MultiField.restrict_to_field(fe,free_dofs,field)
    uhi = FEFunction(U, free_dofs_i, zeros(num_dirichlet_dofs(U)))
    push!(blocks,uhi)
  end
  Gridap.MultiField.MultiFieldFEFunction(free_dofs,fe,blocks)
end
