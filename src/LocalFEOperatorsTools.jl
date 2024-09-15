

function _to_trial_basis(v)
  v
end

function _to_trial_basis(v::Gridap.FESpaces.SingleFieldFEBasis)
  if  Gridap.FESpaces.BasisStyle(v) == Gridap.FESpaces.TestBasis()
    cell_basis=lazy_map(transpose,v.cell_basis)
    return Gridap.FESpaces.SingleFieldFEBasis(cell_basis,
                                              v.trian,
                                              Gridap.FESpaces.TrialBasis(),
                                              v.domain_style)
  else
    return v
  end
end

function _to_trial_basis(v::NTuple{N,Gridap.FESpaces.SingleFieldFEBasis}) where N
  map(_to_trial_basis,v)
end

function _to_trial_basis(v::Gridap.MultiField.MultiFieldCellField)
  if _get_basis_style(v) == Gridap.FESpaces.TestBasis()
    nfields = length(v.single_fields)
    all_febases = Gridap.MultiField.MultiFieldFEBasisComponent[]
    for field_i in 1:nfields
      dv_i = v.single_fields[field_i].single_field
      @assert Gridap.FESpaces.BasisStyle(dv_i) == Gridap.FESpaces.TestBasis()
      dv_i = _to_trial_basis(dv_i)
      dv_i_b = Gridap.MultiField.MultiFieldFEBasisComponent(dv_i,field_i,nfields)
      push!(all_febases,dv_i_b)
    end
    Gridap.MultiField.MultiFieldCellField(all_febases)
  else
    return v
  end
end

function _evaluate_forms(op,v)
  lhs_skeleton=GridapHybrid._find_skeleton(op.LHS_contribs)
  if length(lhs_skeleton)==1
    LHS_contribs = GridapHybrid._merge_bulk_and_skeleton_matrix_contributions(op.LHS_contribs)
  else
    Gridap.Helpers.@check length(lhs_skeleton)==0
    LHS_contribs = op.LHS_contribs
  end
  vs = _get_test_fe_basis(op)
  RHS_contribs = op.RHS_form(_to_trial_basis(v),vs)
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

# Basis of e.g. MultiFieldFESpace
function _get_basis_style(v::Gridap.FESpaces.FEBasis)
  Gridap.FESpaces.BasisStyle(v)
end

function _get_basis_style(v)
  nothing
end


function _get_basis_style(v::NTuple{N,<:CellField}) where N
  bs=map(_get_basis_style,v)
  bs[findfirst(x->x!=nothing,bs)]
end

function _get_basis_style(v::Gridap.CellData.OperationCellField)
  bs=map(_get_basis_style,v.args)
  first=findfirst(x->x!=nothing,bs)
  first!= nothing ? bs[first] : nothing
end

# Basis of a MultiFieldFESpace
function _get_basis_style(v::Gridap.MultiField.MultiFieldCellField)
  Gridap.FESpaces.BasisStyle(first(v))
end

function _get_basis_style(v::NTuple{N,<:Gridap.FESpaces.FEBasis}) where N
  first(map(Gridap.FESpaces.BasisStyle,v))
end

function _convert_cell_dofs_to_facet_dofs(O,cell_dofs)
  t  = get_triangulation(O)
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
                 get_cell_dof_ids(O))[1]
  end
  cell_dofs
end

function _is_on_skeleton(v::Gridap.CellData.CellField)
  t  = get_triangulation(v)
  m  = get_background_model(t)
  Dt = num_cell_dims(t)
  Dm = num_cell_dims(m)
  if (Dt!=Dm)
    @notimplementedif (Dt+1) != Dm
    return true
  else
    return false
  end
end

function _num_funs_x_cell(v::Gridap.FESpaces.SingleFieldFEBasis)
  if (_is_on_skeleton(v))
    m  = get_background_model(get_triangulation(v))
    cell_wise_facets = _get_cell_wise_facets(m)
    n = 0
    cfs=cell_wise_facets[1]
    for f in cfs
      n = n + length(v.cell_basis[f])
    end
    return n
  else
    return length(v.cell_basis[1])
  end
end

function _cell_vector_facets_split(v::Gridap.FESpaces.SingleFieldFEBasis)
  Gridap.Helpers.@check _is_on_skeleton(v)
  result=UnitRange{Int64}[]
  m  = get_background_model(get_triangulation(v))
  cell_wise_facets = _get_cell_wise_facets(m)
  s   = 1
  cfs = cell_wise_facets[1]
  for f in cfs
    n = length(v.cell_basis[f])
    e = s+n-1
    push!(result,s:e)
    s = s + n
  end
  return result
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

function _generate_image_space_span(op,
                                    O::Gridap.FESpaces.SingleFieldFESpace,
                                    v::Gridap.MultiField.MultiFieldCellField,
                                    cell_dofs,
                                    basis_style)
  all_febases = Gridap.MultiField.MultiFieldFEBasisComponent[]
  nfields     = length(v.single_fields)
  s           = 1
  O_basis       = get_fe_basis(O)
  O_basis_data  = Gridap.CellData.get_data(O_basis)
  trian_O_basis = get_triangulation(O_basis)
  for i=1:nfields
    field=v.single_fields[i].single_field
    n = _num_funs_x_cell(field)
    e = s + n - 1
    view_range=s:e
    cell_dofs_current_field=lazy_map(x->view(x,:,view_range),cell_dofs)
    cell_data = lazy_map(Gridap.Fields.linear_combination,
                         cell_dofs_current_field,
                         O_basis_data)
    if (basis_style==Gridap.FESpaces.TrialBasis())
      cell_data = lazy_map(transpose,cell_data)
    end
    du_i=Gridap.FESpaces.similar_fe_basis(O_basis,
                                     cell_data,
                                     trian_O_basis,
                                     basis_style,
                                     DomainStyle(O_basis))
    du_i_b = Gridap.MultiField.MultiFieldFEBasisComponent(du_i,i,nfields)
    push!(all_febases,du_i_b)
    s = e+1
  end
  Gridap.MultiField.MultiFieldCellField(all_febases)
end

function _generate_image_space_span(O::Gridap.FESpaces.SingleFieldFESpace,
                                    cell_dofs,
                                    basis_style)
  cell_dofs=_convert_cell_dofs_to_facet_dofs(O,cell_dofs)
  basis = get_fe_basis(O)
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