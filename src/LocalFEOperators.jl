# TO-DO: better names?
# TO-DO: different implementation strategy? Traits versus ...

abstract type FieldTypeAtCommonFaces end;
struct SingleValued <: FieldTypeAtCommonFaces end;
struct MultiValued  <: FieldTypeAtCommonFaces end;

struct LocalFEOperator{FT<:FieldTypeAtCommonFaces,T1,T2,T3,T4,T5,T6} <: GridapType
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
                           test_space_ds_decomp  :: Union{Nothing,<:MultiFieldFESpace}=nothing,
                           field_type_at_common_faces::FieldTypeAtCommonFaces=SingleValued()) where {A,B}

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

     FT=typeof(field_type_at_common_faces)
     T1=typeof(LHS_form)
     T2=typeof(RHS_form)
     T3=typeof(trial_space)
     T4=typeof(test_space)
     T5=typeof(trial_space_ds_decomp)
     T6=typeof(test_space_ds_decomp)
     new{FT,T1,T2,T3,T4,T5,T6}(LHS_form,
                      RHS_form,
                      trial_space,
                      test_space,
                      trial_space_ds_decomp,
                      test_space_ds_decomp,
                      LHS_contribs)
  end
end

function _get_test_fe_basis(op::LocalFEOperator{FT,A,B,C,D,Nothing}) where {FT,A,B,C,D}
  get_fe_basis(op.test_space)
end

function _get_test_fe_basis(op::LocalFEOperator{FT,A,B,C,D,<:MultiFieldFESpace}) where {FT,A,B,C,D}
  get_fe_basis(op.test_space_ds_decomp)
end

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

function _evaluate_forms(op::LocalFEOperator,v)
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

function _generate_image_space_span(op::LocalFEOperator,
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

function _generate_image_space_span(op::LocalFEOperator{<:SingleValued},
                                    O::Gridap.MultiField.MultiFieldFESpace,
                                    v::Gridap.CellData.CellField,
                                    cell_dofs,
                                    basis_style) where N
  cell_dofs_field_offsets=_compute_cell_dofs_field_offsets(O)
  all_febases = Gridap.MultiField.MultiFieldFEBasisComponent[]
  nfields = length(O.spaces)
  for i=1:nfields
    view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
    cell_dofs_current_field=lazy_map(x->view(x,view_range,:),cell_dofs)
    du_i=_generate_image_space_span(O.spaces[i],cell_dofs_current_field,basis_style)
    du_i_b = Gridap.MultiField.MultiFieldFEBasisComponent(du_i,i,nfields)
    push!(all_febases,du_i_b)
  end
  Gridap.MultiField.MultiFieldCellField(all_febases)
end

function _generate_image_space_span(op::LocalFEOperator{<:MultiValued},
                                    O::Gridap.MultiField.MultiFieldFESpace,
                                    v::Gridap.MultiField.MultiFieldCellField,
                                    cell_dofs,
                                    basis_style) where N

  nfields = length(O.spaces)
  if (basis_style==Gridap.FESpaces.TrialBasis())
    mf_basis=get_trial_fe_basis(O)
  else
    mf_basis=get_fe_basis(O)
  end
  multi_fields = Gridap.MultiField.MultiFieldCellField[]
  
  col_stride=0
  for i=1:nfields
    field=v.single_fields[i].single_field
    n = _num_funs_x_cell(field)
    col_stride+=n 
  end 
  
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
      
      sj=(i-1)*col_stride+1
      for j=1:nfields
        field=v.single_fields[j].single_field
        nj = _num_funs_x_cell(field)
        ej = sj + nj - 1

        if (_is_on_skeleton(field))
          stride=sum([length(k) for k in _cell_vector_facets_split(field)])
          col_dofs_split=UnitRange{Int64}[]
          cs = sj 
          for k=1:length(_cell_vector_facets_split(field))
            ce = cs+stride-1
            push!(col_dofs_split,cs:ce)
            cs=cs+col_stride
          end
        else 
          col_dofs_split=sj:ej
        end 
        
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

function (op::LocalFEOperator)(v::Union{Gridap.MultiField.MultiFieldCellField,
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
function (op::LocalFEOperator)(v::FEFunction)
  LHSf,RHSf = _evaluate_forms(op,v)
  cell_dofs = _compute_cell_dofs(LHSf,RHSf)
  _generate_fe_function(op,cell_dofs)
end
