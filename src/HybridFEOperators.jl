
mutable struct HybridFEMatrix{T,TB,TS} <: AbstractMatrix{T}
  bulk_fields        :: TB
  skeleton_fields    :: TS
  trial_hybridizable :: MultiFieldFESpace
  test_hybridizable  :: MultiFieldFESpace
  trial_skeleton     :: FESpace
  test_skeleton      :: FESpace
  mat_contribs       :: DomainContribution # Cell-wise un-assembled matrix contributions
end

mutable struct HybridFEVector{T,V,TB,TS} <: AbstractVector{T}
  bulk_fields        :: TB
  skeleton_fields    :: TS
  trial_hybridizable :: MultiFieldFESpace
  test_hybridizable  :: MultiFieldFESpace
  trial_skeleton     :: FESpace
  test_skeleton      :: FESpace
  vec_contribs       :: DomainContribution # Cell-wise un-assembled vector contributions
  vec                :: V
  function HybridFEVector(bulk_fields::TB,
                          skeleton_fields::TS,
                          trial_hybridizable::MultiFieldFESpace,
                          test_hybridizable::MultiFieldFESpace,
                          trial_skeleton::FESpace,
                          test_skeleton::FESpace,
                          vec_contribs::DomainContribution) where {TB,TS}
    vec=assemble_vector(vec_contribs,test_hybridizable)
    T=Float64
    V=typeof(vec)
    new{T,V,TB,TS}(bulk_fields,
                   skeleton_fields,
                   trial_hybridizable,
                   test_hybridizable,
                   trial_skeleton,
                   test_skeleton,
                   vec_contribs,
                   vec)
  end
end

function Base.copy(m::HybridFEMatrix)
  # TO-DO: this is not a deep copy
  m
end

function Base.size(v::HybridFEMatrix)
  (num_free_dofs(v.test_hybridizable),num_free_dofs(v.trial_hybridizable))
end

function Base.copy(v::HybridFEVector)
  # TO-DO: this is not a deep copy
  v
end

function Base.size(v::HybridFEVector)
  (num_free_dofs(v.test_hybridizable),)
end

function Base.getindex(v::HybridFEVector,i)
  v.vec[i]
end

function _update_vec!(v::HybridFEVector)
  Gridap.FESpaces.assemble_vector!(v.vec_contribs,v.vec,v.test_hybridizable)
end

struct HybridFEOperator{TB,TS} <: FEOperator
  res::Function
  jac::Function
  trial::MultiFieldFESpace
  test::MultiFieldFESpace
  bulk_fields::TB
  skeleton_fields::TS
  trial_skeleton::FESpace
  test_skeleton::FESpace
  assem_hybridizable::SparseMatrixAssembler
  assem_skeleton::SparseMatrixAssembler
  function HybridFEOperator(
    res::Function,
    jac::Function,
    trial::MultiFieldFESpace,
    test::MultiFieldFESpace,
    bulk_fields::TB,
    skeleton_fields::TS,
    trial_skeleton::FESpace,
    test_skeleton::FESpace,
    assem_hybridizable::SparseMatrixAssembler,
    assem_skeleton::SparseMatrixAssembler) where {TB,TS}
    @assert ! isa(test,TrialFESpace) """\n
    It is not allowed to build a FEOperator with a test space of type TrialFESpace.

    Make sure that you are writing first the trial space and then the test space when
    building a HybridFEOperator.
    """
    new{TB,TS}(res,
               jac,
               trial,
               test,
               bulk_fields,
               skeleton_fields,
               trial_skeleton,
               test_skeleton,
               assem_hybridizable,
               assem_skeleton)
  end
end

function HybridFEOperator(
  res::Function,
  jac::Function,
  trial::MultiFieldFESpace,
  test::MultiFieldFESpace,
  bulk_fields::Vector{<:Integer},
  skeleton_fields::Vector{<:Integer})
  assem_hybridizable=SparseMatrixAssembler(trial,test)
  trial_skeleton,test_skeleton=_setup_fe_spaces_skeleton_system(trial,
                                                                   test,
                                                                   skeleton_fields)
  assem_skeleton=SparseMatrixAssembler(trial_skeleton,test_skeleton)
  HybridFEOperator(res,
                    jac,
                    trial,
                    test,
                    bulk_fields,
                    skeleton_fields,
                    trial_skeleton,
                    test_skeleton,
                    assem_skeleton,
                    assem_hybridizable)
end

function HybridFEOperator(
  res::Function,
  trial::MultiFieldFESpace,
  test::MultiFieldFESpace,
  bulk_fields::Vector{<:Integer},
  skeleton_fields::Vector{<:Integer})
  jac(u,du,dv) = jacobian(x->res(x,dv),u)
  HybridFEOperator(res,jac,trial,test,bulk_fields,skeleton_fields)
end

Gridap.FESpaces.get_test(op::HybridFEOperator) = op.test
Gridap.FESpaces.get_trial(op::HybridFEOperator) = op.trial

function Gridap.FESpaces.allocate_residual(op::HybridFEOperator,uh)
  V = Gridap.FESpaces.get_test(op)
  v = Gridap.FESpaces.get_fe_basis(V)
  dcres=op.res(uh,v)
  TB=typeof(op.bulk_fields)
  TS=typeof(op.skeleton_fields)
  HybridFEVector(op.bulk_fields,
                 op.skeleton_fields,
                 op.trial,
                 op.test,
                 op.trial_skeleton,
                 op.test_skeleton,
                 dcres)
end

function Gridap.FESpaces.residual!(b::HybridFEVector,op::HybridFEOperator,uh)
  V = Gridap.FESpaces.get_test(op)
  v = Gridap.FESpaces.get_fe_basis(V)
  dcres=op.res(uh,v)
  b.vec_contribs=dcres
  _update_vec!(b)
  b
end

function Gridap.FESpaces.allocate_jacobian(op::HybridFEOperator,uh)
  U = Gridap.FESpaces.get_trial(op)
  V = Gridap.FESpaces.get_test(op)
  du = Gridap.FESpaces.get_trial_fe_basis(U)
  v = Gridap.FESpaces.get_fe_basis(V)
  dcjac=op.jac(uh,du,v)
  TB=typeof(op.bulk_fields)
  TS=typeof(op.skeleton_fields)
  HybridFEMatrix{Float64,TB,TS}(op.bulk_fields,
                                op.skeleton_fields,
                                op.trial,
                                op.test,
                                op.trial_skeleton,
                                op.test_skeleton,
                                dcjac)
end

function Gridap.FESpaces.jacobian!(A::HybridFEMatrix,op::HybridFEOperator,uh)
  U = Gridap.FESpaces.get_trial(op)
  V = Gridap.FESpaces.get_test(op)
  du = Gridap.FESpaces.Gridap.FESpaces.get_trial_fe_basis(U)
  v = Gridap.FESpaces.get_fe_basis(V)
  dcjac=op.jac(uh,du,v)
  A.mat_contribs=dcjac
  A
end

function Gridap.FESpaces.residual_and_jacobian!(
  b::HybridFEVector,A::HybridFEMatrix,op::HybridFEOperator,uh)

  Gridap.Helpers.@check b.test_skeleton === op.test_skeleton
  Gridap.Helpers.@check A.trial_skeleton === op.trial_skeleton

  U = Gridap.FESpaces.get_trial(op)
  V = Gridap.FESpaces.get_test(op)
  du = Gridap.FESpaces.get_trial_fe_basis(U)
  v = Gridap.FESpaces.get_fe_basis(V)
  dcjac=op.jac(uh,du,v)
  dcres=op.res(uh,v)
  b.vec_contribs=dcres
  _update_vec!(b)
  A.mat_contribs=dcjac
  (b,A)
end

function Gridap.FESpaces.residual_and_jacobian(op::HybridFEOperator,uh)
  U = Gridap.FESpaces.get_trial(op)
  V = Gridap.FESpaces.get_test(op)
  du = Gridap.FESpaces.get_trial_fe_basis(U)
  v = Gridap.FESpaces.get_fe_basis(V)
  dcjac=op.jac(uh,du,v)
  dcres=op.res(uh,v)
  TB=typeof(op.bulk_fields)
  TS=typeof(op.skeleton_fields)
  A=HybridFEMatrix{Float64,TB,TS}(op.bulk_fields,
                                  op.skeleton_fields,
                                  op.trial,
                                  op.test,
                                  op.trial_skeleton,
                                  op.test_skeleton,
                                  dcjac)
  b=HybridFEVector(op.bulk_fields,
                   op.skeleton_fields,
                   op.trial,
                   op.test,
                   op.trial_skeleton,
                   op.test_skeleton,
                   dcres)
  (b, A)
end
