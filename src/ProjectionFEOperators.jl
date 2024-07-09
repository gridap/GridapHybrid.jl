
struct ProjectionFEOperator{T1,T2,T3,T4} <: GridapType
  LHS_form               :: T1
  RHS_form               :: T2
  trial_space            :: T3
  test_space             :: T4
  LHS_contribs           :: Gridap.CellData.DomainContribution

  function ProjectionFEOperator(weakform,
                           trial_space :: A,
                           test_space  :: B) where {A,B}

     us = get_trial_fe_basis(trial_space)
     vs = get_fe_basis(test_space)

     LHS_form, RHS_form = weakform
     LHS_contribs       = LHS_form(us,vs)

     T1=typeof(LHS_form)
     T2=typeof(RHS_form)
     T3=typeof(trial_space)
     T4=typeof(test_space)
     new{T1,T2,T3,T4}(LHS_form,
                      RHS_form,
                      trial_space,
                      test_space,
                      LHS_contribs)
  end
end

function _get_test_fe_basis(op::ProjectionFEOperator)
  get_fe_basis(op.test_space)
end

function _generate_image_space_span(op::ProjectionFEOperator,
                                    O::Gridap.MultiField.MultiFieldFESpace,
                                    v::Gridap.MultiField.MultiFieldCellField,
                                    cell_dofs,
                                    basis_style)

  nfields = length(O.spaces)
  if (basis_style==Gridap.FESpaces.TrialBasis())
    mf_basis=get_trial_fe_basis(O)
  else
    mf_basis=get_fe_basis(O)
  end
  multi_fields = Gridap.MultiField.MultiFieldCellField[]
  
  s=1
  for i=1:nfields
    field=v.single_fields[i].single_field
    n = _num_funs_x_cell(field)
    e = s + n - 1
    row_view_range=s:e
    if (_is_on_skeleton(field))
      single_fields = Gridap.CellData.GenericCellField[]
      m = get_background_model(get_triangulation(field))
      cell_wise_facets_ids = _get_cell_wise_facets(m)
      tskel = Skeleton(m)


      row_dofs_split=[ (s-1).+i for i in _cell_vector_facets_split(field)]
      sj=1
      for j=1:nfields
        field=v.single_fields[j].single_field
        nj = _num_funs_x_cell(field)
        ej = sj + nj - 1
        col_dofs_split=sj:ej
        skel_facet_dofs=SkeletonVectorFromSplitDoFsCellVector(tskel.glue,
                                              cell_wise_facets_ids,
                                              cell_dofs,
                                              row_dofs_split,
                                              col_dofs_split)

        facet_shapefuns =
           Gridap.MultiField.MultiFieldFEBasisComponent(mf_basis[i].single_field,j,nfields)
        facet_shapefuns=Gridap.CellData.get_data(facet_shapefuns)

        skel_facet_shapefuns=_transform_face_to_cell_lface_array(tskel.glue,facet_shapefuns;
                                                          add_naive_innermost_block_level=true)
        cell_field=lazy_map(linear_combination,skel_facet_dofs,skel_facet_shapefuns)
        cf = Gridap.CellData.GenericCellField(cell_field, tskel, DomainStyle(mf_basis[i]))
        push!(single_fields,cf)
        sj=ej+1
      end
      multi_field=Gridap.MultiField.MultiFieldCellField(single_fields)
    else
       single_fields = Gridap.MultiField.MultiFieldFEBasisComponent[]
       sj=1
       for j=1:nfields
          field=v.single_fields[j].single_field
          nj = _num_funs_x_cell(field)
          ej = sj + nj - 1
          col_view_range=sj:ej
          cell_dofs_current_field=lazy_map(x->view(x,row_view_range,col_view_range),cell_dofs)
          du_j = _generate_image_space_span(O.spaces[i],cell_dofs_current_field,basis_style)
          du_j_b = Gridap.MultiField.MultiFieldFEBasisComponent(du_j,j,nfields)
          push!(single_fields,du_j_b)
          sj=ej+1
       end
       multi_field=Gridap.MultiField.MultiFieldCellField(single_fields)
    end
    s = e+1
    push!(multi_fields,multi_field)
  end
  multi_fields
end

function _generate_cell_field(field_id, nfields, U::Gridap.FESpaces.SingleFieldFESpace,cell_dofs)
  t  = get_triangulation(U)
  m  = get_background_model(t)
  Dt = num_cell_dims(t)
  Dm = num_cell_dims(m)
  skel=Skeleton(m)
  if (Dt!=Dm)
    @notimplementedif (Dt+1) != Dm
    # facet_dofs=_convert_cell_dofs_to_facet_dofs(U,cell_dofs)
    # skel_facet_dofs=SkeletonVectorFromFacetVector(skel.glue,cell_wise_facets_ids,facet_dofs)

    cell_wise_facets_ids = _get_cell_wise_facets(m)

    basis=get_fe_basis(U)
    facet_shapefuns=Gridap.CellData.get_data(basis)
    row_dofs_split=_cell_vector_facets_split(basis)
    col_dofs_split=1

    skel_facet_dofs=SkeletonVectorFromSplitDoFsCellVector(skel.glue,
                                      cell_wise_facets_ids,
                                      cell_dofs,
                                      row_dofs_split,
                                      col_dofs_split)

    skel_facet_shapefuns=_transform_face_to_cell_lface_array(skel.glue,facet_shapefuns;
                                                             add_naive_innermost_block_level=true)

    cell_field=lazy_map(linear_combination,skel_facet_dofs,skel_facet_shapefuns)
    cf = Gridap.CellData.GenericCellField(cell_field, skel, DomainStyle(get_fe_basis(U)))
  else 
    free_dofs = Gridap.FESpaces.gather_free_values(U,cell_dofs)
    cell_field = Gridap.CellData.get_data(FEFunction(U,free_dofs))
    cf = Gridap.CellData.GenericCellField(cell_field, t, ReferenceDomain())
  end
  cf
end

function _generate_cell_field(op::ProjectionFEOperator,cell_dofs)
  U = op.trial_space
  cell_dofs_field_offsets=_compute_cell_dofs_field_offsets(U)
  single_fields = Gridap.CellData.GenericCellField[]
  nfields = length(U.spaces)
  for i=1:nfields
      Ui=U.spaces[i]
      view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
      cell_dofs_current_field=lazy_map(x->view(x,view_range),cell_dofs)
      du_i=_generate_cell_field(i, nfields, U.spaces[i],cell_dofs_current_field)
      push!(single_fields,du_i)
   end
   Gridap.MultiField.MultiFieldCellField(single_fields)
end

function (op::ProjectionFEOperator)(v::Union{Gridap.MultiField.MultiFieldCellField,
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
function (op::ProjectionFEOperator)(v::FEFunction)
  LHSf,RHSf = _evaluate_forms(op,v)
  cell_dofs = _compute_cell_dofs(LHSf,RHSf)
  _generate_cell_field(op,cell_dofs)
end
