struct LocalFEOperator{T1,T2,T3,T4,T5,T6} <: GridapType
  LHS_form               :: T1
  RHS_form               :: T2
  trial_space            :: T3
  test_space             :: T4
  trial_space_ds_decomp  :: T5 # Direct sum decomposition of trial space (optional)
  test_space_ds_decomp   :: T6 # Direct sum decomposition of test space  (optional)
  LHS_contribs           :: Gridap.CellData.DomainContribution

  function LocalFEOperator(weakform,
                           trial_space :: A,
                           test_space  :: B;
                           trial_space_ds_decomp :: Union{Nothing,<:MultiFieldFESpace}=nothing,
                           test_space_ds_decomp  :: Union{Nothing,<:MultiFieldFESpace}=nothing) where {A,B}

    if trial_space_ds_decomp != nothing
      Gridap.Helpers.@check test_space_ds_decomp != nothing
      Gridap.Helpers.@check length(trial_space_ds_decomp) == length(trial_space_ds_decomp)
    end

    if (trial_space_ds_decomp == nothing)
        us = get_trial_fe_basis(trial_space)
        vs = get_fe_basis(test_space)
    else
        us = get_trial_fe_basis(trial_space_ds_decomp)
        vs = get_fe_basis(test_space_ds_decomp)
    end

     LHS_form, RHS_form = weakform
     LHS_contribs       = LHS_form(us,vs)

     T1=typeof(LHS_form)
     T2=typeof(RHS_form)
     T3=typeof(trial_space)
     T4=typeof(test_space)
     T5=typeof(trial_space_ds_decomp)
     T6=typeof(test_space_ds_decomp)
     new{T1,T2,T3,T4,T5,T6}(LHS_form,
                      RHS_form,
                      trial_space,
                      test_space,
                      trial_space_ds_decomp,
                      test_space_ds_decomp,
                      LHS_contribs)
  end
end

function _get_test_fe_basis(op::LocalFEOperator{A,B,C,D,Nothing}) where {A,B,C,D}
  get_fe_basis(op.test_space)
end

function _get_test_fe_basis(op::LocalFEOperator{A,B,C,D,<:MultiFieldFESpace}) where {A,B,C,D}
  get_fe_basis(op.test_space_ds_decomp)
end

function _evaluate_forms(op::LocalFEOperator,v)
  lhs_skeleton=GridapHybrid._find_skeleton(op.LHS_contribs)
  if length(lhs_skeleton)==1
    LHS_contribs = Gridap.Hybrid._merge_bulk_and_skeleton_contributions(op.LHS_contribs)
  else
    Gridap.Helpers.@check length(lhs_skeleton)==0
    LHS_contribs = op.LHS_contribs
  end
  vs = _get_test_fe_basis(op)
  RHS_contribs = op.RHS_form(v,vs)
  rhs_skeleton=GridapHybrid._find_skeleton(RHS_contribs)
  if length(rhs_skeleton)==1
    RHS_contribs = GridapHybrid._merge_bulk_and_skeleton_matrix_contributions(RHS_contribs)
  end
  LHS_contribs, RHS_contribs
end

function _compute_cell_dofs(LHS_contribs,RHS_contribs)
  A_array = lazy_map(DensifyInnerMostBlockLevelMap(),get_array(LHS_contribs))
  B_array = lazy_map(DensifyInnerMostBlockLevelMap(),get_array(RHS_contribs))
  dofs = lazy_map(\,A_array,B_array)
end


function (op::LocalFEOperator)(v::Gridap.MultiField.MultiFieldCellField)
  basis_style = Gridap.FESpaces.BasisStyle(first(v))

  Gridap.Helpers.@notimplementedif basis_style == Gridap.FESpaces.TestBasis()

  LHSf,RHSf= _evaluate_forms(op,v)
  cell_dofs= _compute_cell_dofs(LHSf,RHSf)

  # We get the test space so that we don't deal with
  # the transpose underlying the trial basis
  test_basis = get_fe_basis(op.test_space)
  cell_data =  lazy_map(Gridap.Fields.linear_combination,
                        cell_dofs,
                        Gridap.CellData.get_data(test_basis))

  if (basis_style==Gridap.FESpaces.TrialBasis())
    cell_data = lazy_map(transpose,cell_data)
  end

  Gridap.FESpaces.similar_fe_basis(test_basis,
                                   cell_data,
                                   Gridap.CellData.get_domains(LHSf)...,
                                   basis_style,
                                   DomainStyle(test_basis))
end

function (op::LocalFEOperator)(v::Gridap.MultiField.MultiFieldFEFunction)
  LHSf,RHSf = _evaluate_forms(op,v)
  cell_dofs = _compute_cell_dofs(LHSf,RHSf)
  U         = op.trial_space
  free_dofs = Gridap.FESpaces.gather_free_values(U,cell_dofs)
  FEFunction(U,free_dofs)
end
