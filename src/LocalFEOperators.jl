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
    LHS_contribs = GridapHybrid._merge_bulk_and_skeleton_matrix_contributions(op.LHS_contribs)
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

function _get_basis_style(v::Gridap.FESpaces.FEBasis)
  Gridap.FESpaces.BasisStyle(v)
end

function _get_basis_style(v::Gridap.MultiField.MultiFieldCellField)
  Gridap.FESpaces.BasisStyle(first(v))
end

function _generate_image_space_span(U::Gridap.FESpaces.SingleFieldFESpace,cell_dofs,basis_style)
  t  = get_triangulation(U)
  m  = get_background_model(t)
  Dt = num_cell_dims(t)
  Dm = num_cell_dims(m)
  if (Dt!=Dm)
    @notimplementedif (Dt+1) != Dm
    cell_wise_facets = _get_cell_wise_facets(m)
    cells_around_facets = _get_cells_around_facets(m)
    cell_dofs=convert_cell_wise_dofs_array_to_facet_dofs_array(
                 cells_around_facets,
                 cell_wise_facets,
                 cell_dofs,
                 get_cell_dof_ids(U))[1]
  end
  basis = get_fe_basis(U)
  cell_data = lazy_map(Gridap.Fields.linear_combination,
                       cell_dofs,
                       Gridap.CellData.get_data(basis))
  if (basis_style==Gridap.FESpaces.TrialBasis())
    cell_data = lazy_map(transpose,cell_data)
  end
  Gridap.FESpaces.similar_fe_basis(basis,
                                  cell_data,
                                  get_triangulation(basis),
                                  basis_style,
                                  DomainStyle(basis))
end

function _compute_cell_dofs_field_offsets(U::Gridap.MultiField.MultiFieldFESpace)
  nfields = length(U.spaces)
  cell_dofs_field_offsets=Vector{Int}(undef,nfields+1)
  cell_dofs_field_offsets[1]=1
  for i=1:nfields
    Ui = U.spaces[i]
    t = get_triangulation(Ui)
    m = get_background_model(t)
    Dt,Dm = num_cell_dims(t),num_cell_dims(m)
    if (Dt!=Dm)
     Ui_dofs = _get_facet_dofs_cell_wise(Ui)
    else
     Ui_dofs = get_cell_dof_ids(Ui)
    end
    cell_dofs_field_offsets[i+1]=cell_dofs_field_offsets[i]+length(Ui_dofs[1])
  end
  cell_dofs_field_offsets
end

function _generate_image_space_span(U::Gridap.MultiField.MultiFieldFESpace,cell_dofs,basis_style)
  cell_dofs_field_offsets=_compute_cell_dofs_field_offsets(U)
  all_febases = Gridap.MultiField.MultiFieldFEBasisComponent[]
  nfields = length(U.spaces)
  for i=1:nfields
    view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
    cell_dofs_current_field=lazy_map(x->view(x,view_range,:),cell_dofs)
    du_i=_generate_image_space_span(U.spaces[i],cell_dofs_current_field,basis_style)
    du_i_b = Gridap.MultiField.MultiFieldFEBasisComponent(du_i,i,nfields)
    push!(all_febases,du_i_b)
  end
  Gridap.MultiField.MultiFieldCellField(all_febases)
end

function _generate_image_space_span(op::LocalFEOperator,cell_dofs,basis_style)
  # We get the test space so that we don't deal with
  # the transpose underlying the trial basis
  U = op.test_space
  _generate_image_space_span(U,cell_dofs,basis_style)
end

function _generate_fe_function(U::Gridap.FESpaces.SingleFieldFESpace,cell_dofs)
  t  = get_triangulation(U)
  m  = get_background_model(t)
  Dt = num_cell_dims(t)
  Dm = num_cell_dims(m)
  if (Dt!=Dm)
    @notimplementedif (Dt+1) != Dm
    cell_wise_facets = _get_cell_wise_facets(m)
    cells_around_facets = _get_cells_around_facets(m)
    cell_dofs=convert_cell_wise_dofs_array_to_facet_dofs_array(
                 cells_around_facets,
                 cell_wise_facets,
                 cell_dofs,
                 get_cell_dof_ids(U))[1]
  end
  free_dofs = Gridap.FESpaces.gather_free_values(U,cell_dofs)
  FEFunction(U,free_dofs)
end

function _generate_fe_function(U::Gridap.MultiField.MultiFieldFESpace,cell_dofs)
  cell_dofs_field_offsets=_compute_cell_dofs_field_offsets(U)
  # TO-DO: extract eltype from U, instead of a hard-coded Float64
  free_values = Vector{Float64}(undef,num_free_dofs(U))
  nfields = length(U.spaces)
  s=1
  for i=1:nfields
    Ui=U.spaces[i]
    e = s + (num_free_dofs(Ui)-1)
    view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
    cell_dofs_current_field=lazy_map(x->view(x,view_range),cell_dofs)
    du_i=_generate_fe_function(U.spaces[i],cell_dofs_current_field)
    free_values[s:e] .= get_free_dof_values(du_i)
    s=e+1
  end
  FEFunction(U, free_values)
end

function _generate_fe_function(op::LocalFEOperator,cell_dofs)
  U = op.trial_space
  _generate_fe_function(U,cell_dofs)
end

function (op::LocalFEOperator)(v::Union{Gridap.FESpaces.FEBasis,
                                        Gridap.MultiField.MultiFieldCellField})
  basis_style = _get_basis_style(v)
  Gridap.Helpers.@notimplementedif basis_style == Gridap.FESpaces.TestBasis()
  LHSf,RHSf= _evaluate_forms(op,v)
  cell_dofs= _compute_cell_dofs(LHSf,RHSf)
  _generate_image_space_span(op,cell_dofs,basis_style)
end

# We want to be as general as possible with the type of objects
# that can appear in place of v. In principle, a given v is valid
# if and only if it is legal to be passed as an argument to
# evaluate_forms
function (op::LocalFEOperator)(v)
  LHSf,RHSf = _evaluate_forms(op,v)
  cell_dofs = _compute_cell_dofs(LHSf,RHSf)
  _generate_fe_function(op,cell_dofs)
end
