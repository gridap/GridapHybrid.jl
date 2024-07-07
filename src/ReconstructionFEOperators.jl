struct ReconstructionFEOperator{T1,T2,T3,T4,T5,T6} <: GridapType
  LHS_form               :: T1
  RHS_form               :: T2
  trial_space            :: T3
  test_space             :: T4
  trial_space_ds_decomp  :: T5 # Direct sum decomposition of trial space (optional)
  test_space_ds_decomp   :: T6 # Direct sum decomposition of test space  (optional)
  LHS_contribs           :: Gridap.CellData.DomainContribution

  function ReconstructionFEOperator(weakform,
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

function _get_test_fe_basis(op::ReconstructionFEOperator{FT,A,B,C,D,Nothing}) where {FT,A,B,C,D}
  get_fe_basis(op.test_space)
end

function _get_test_fe_basis(op::ReconstructionFEOperator{FT,A,B,C,D,<:MultiFieldFESpace}) where {FT,A,B,C,D}
  get_fe_basis(op.test_space_ds_decomp)
end

function _generate_image_space_span(op::ReconstructionFEOperator,
                                    image_space::Gridap.MultiField.MultiFieldFESpace,
                                    v::Gridap.MultiField.MultiFieldCellField,
                                    cell_dofs,
                                    basis_style)

  cell_dofs_field_offsets=_compute_cell_dofs_field_offsets(image_space)
  all_febases = Vector{Gridap.MultiField.MultiFieldFEBasisComponent}(undef,2)
  
  row_view_range=cell_dofs_field_offsets[1]:cell_dofs_field_offsets[2]-1
  
  nfields = length(v.single_fields)
  field=v.single_fields[1].single_field
  n=_num_funs_x_cell(field)
  col_view_range_1=1:n 
  cell_dofs_current_field_1=lazy_map(x->view(x,row_view_range,col_view_range_1),cell_dofs)
  du_i_1=_generate_image_space_span(image_space.spaces[1],cell_dofs_current_field_1,basis_style)
  du_i_b_1=Gridap.MultiField.MultiFieldFEBasisComponent(du_i_1,1,nfields)
  all_febases[1]=du_i_b_1
  field=v.single_fields[2].single_field

  # Had to use another local variable for col_view_range to avoid 
  # changing the reference in the closure above. 
  col_view_range_2=n+1:n+_num_funs_x_cell(field)
  cell_dofs_current_field_2=lazy_map(x->view(x,row_view_range,col_view_range_2),cell_dofs)
  du_i_2=_generate_image_space_span(image_space.spaces[1],cell_dofs_current_field_2,basis_style)
  du_i_b_2 = Gridap.MultiField.MultiFieldFEBasisComponent(du_i_2,2,nfields)
  all_febases[2]=du_i_b_2
  Gridap.MultiField.MultiFieldCellField(all_febases)  
end

function _generate_cell_field(op::ReconstructionFEOperator,cell_dofs)
    U=op.trial_space
    all_components = Vector{GenericCellField}(undef,2)
    cell_dofs_field_offsets=_compute_cell_dofs_field_offsets(U)
    view_range=cell_dofs_field_offsets[1]:cell_dofs_field_offsets[2]-1
    cell_dofs_current_field=lazy_map(x->view(x,view_range),cell_dofs)
    free_dofs = Gridap.FESpaces.gather_free_values(U[1],cell_dofs_current_field)
    uh=FEFunction(U[1],free_dofs)

    # Replicate the reconstructed function in both blocks
    # To some extent, in my view, this is quite dirty and may be
    # an indication that the current approach that we chose to implement 
    # HHO might rot.
    data=lazy_map(BlockMap(2,1),Gridap.CellData.get_data(uh))
    cf = Gridap.CellData.GenericCellField(data, 
                                          get_triangulation(U[1]), 
                                          ReferenceDomain())
    all_components[1]=cf

    data=lazy_map(BlockMap(2,2),Gridap.CellData.get_data(uh))
    cf = Gridap.CellData.GenericCellField(data, 
                                          get_triangulation(U[1]), 
                                          ReferenceDomain())
    all_components[2]=cf
    Gridap.MultiField.MultiFieldCellField(all_components)
end

function (op::ReconstructionFEOperator)(v::Union{Gridap.MultiField.MultiFieldCellField,
                                        Gridap.FESpaces.FEBasis,
                                        Tuple{N,<:Gridap.CellData.CellField}}) where N
  basis_style = _get_basis_style(v)
  LHSf,RHSf= _evaluate_forms(op,v)
  cell_dofs= _compute_cell_dofs(LHSf,RHSf)
  # We get the test space so that we don't deal with
  # the transpose underlying the trial basis
  O = op.test_space
  _generate_image_space_span(op,O,v,cell_dofs,basis_style)
end

# We want to be as general as possible with the type of objects
# that can appear in place of v. In principle, a given v is valid
# if and only if it is legal to be passed as an argument to
# evaluate_forms
function (op::ReconstructionFEOperator)(v::FEFunction)
  LHSf,RHSf = _evaluate_forms(op,v)
  cell_dofs = _compute_cell_dofs(LHSf,RHSf)
  _generate_cell_field(op,cell_dofs)
end
