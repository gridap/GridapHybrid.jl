struct HybridAffineFEOperator{TB,TS} <: FEOperator
  weakform::Function
  trial::MultiFieldFESpace
  test::MultiFieldFESpace
  bulk_fields::TB
  skeleton_fields::TS
  condensed_op::AffineFEOperator
end

function HybridAffineFEOperator(
  weakform::Function,
  trial :: MultiFieldFESpace,
  test :: MultiFieldFESpace,
  bulk_fields :: TB,
  skeleton_fields :: TS) where {TB<:Vector{<:Integer},TS<:Vector{<:Integer}}

  # Invoke weak form of the hybridizable system
  u = get_trial_fe_basis(trial)
  v = get_fe_basis(test)
  biform, liform  = weakform(u,v)

  # Transform DomainContribution objects of the hybridizable system into a
  # suitable form for assembling the linear system defined on the skeleton
  # (i.e., the hybrid system)
  obiform, oliform = _hybridrizable_to_hybrid_contributions(biform,liform)

  # Pair LHS and RHS terms associated to TempSkeletonTriangulation
  matvec,mat,vec=Gridap.FESpaces._pair_contribution_when_possible(obiform,oliform)

  # Add StaticCondensationMap to matvec terms
  matvec=_add_static_condensation(matvec,bulk_fields,skeleton_fields)

  if (length(skeleton_fields)==1)
    # Single-field skeleton system
    M = trial[skeleton_fields[1]]
    L = test[skeleton_fields[1]]
  else
    # Multi-field skeleton system
    M = MultiFieldFESpace([trial[i] for i in skeleton_fields])
    L = MultiFieldFESpace([test[i] for i in skeleton_fields])
    matvec=_block_skeleton_system_contributions(matvec,L)
  end
  assem=SparseMatrixAssembler(M,L)

  # Attach strong imposition of Dirichlet BCs
  uhd = zero(M)
  matvec,mat=Gridap.FESpaces._attach_dirichlet(matvec,mat,uhd)

  data = Gridap.FESpaces._collect_cell_matrix_and_vector(M,L,matvec,mat,vec)

  A,b = assemble_matrix_and_vector(assem,data)

  condensed_op=AffineFEOperator(M,L,A,b)
  HybridAffineFEOperator(weakform,trial,test,bulk_fields,skeleton_fields,condensed_op)
end

Gridap.FESpaces.get_test(feop::HybridAffineFEOperator) = feop.test
Gridap.FESpaces.get_trial(feop::HybridAffineFEOperator) = feop.trial
function Gridap.FESpaces.solve(op::HybridAffineFEOperator)
  solver = LinearFESolver()
  solve(solver,op)
end

function Gridap.FESpaces.solve!(uh,solver::LinearFESolver,op::HybridAffineFEOperator, cache)
  # Solve linear system defined on the skeleton
  lh = solve(op.condensed_op)

  # Invoke weak form of the hybridizable system
  u = get_trial_fe_basis(op.trial)
  v = get_fe_basis(op.test)
  biform, liform  = op.weakform(u,v)

  # Transform DomainContribution objects of the hybridizable system into a
  # suitable form for assembling the linear system defined on the skeleton
  # (i.e., the hybrid system)
  obiform, oliform = _hybridrizable_to_hybrid_contributions(biform,liform)

  # Pair LHS and RHS terms associated to TempSkeletonTriangulation
  matvec,_,_=Gridap.FESpaces._pair_contribution_when_possible(obiform,oliform)

  # Convert dof-wise dof values of lh into cell-wise dof values lhₖ
  Γ=first(keys(matvec.dict))
  Gridap.Helpers.@check isa(Γ,TempSkeletonTriangulation)
  lhₖ=get_cell_dof_values(lh,Γ)

  # Compute cell-wise dof values of bulk fields out of lhₖ
  t = matvec.dict[Γ]
  m=BackwardStaticCondensationMap(op.bulk_fields,op.skeleton_fields)
  uhphlhₖ=lazy_map(m,t,lhₖ)

  model=get_background_model(Γ)
  cell_wise_facets=_get_cell_wise_facets(model)
  cells_around_facets=_get_cells_around_facets(model)

  nfields=length(op.bulk_fields)+length(op.skeleton_fields)
  m=Gridap.Fields.BlockMap(nfields,op.skeleton_fields)
  L=Gridap.FESpaces.get_fe_space(lh)
  lhₑ=lazy_map(m,
                 convert_cell_wise_dofs_array_to_facet_dofs_array(
                 cells_around_facets,
                 cell_wise_facets,
                 lhₖ,
                 get_cell_dof_ids(L))...)

  assem = SparseMatrixAssembler(op.trial,op.test)
  lhₑ_dofs=get_cell_dof_ids(op.trial,get_triangulation(L))

  Ω = op.trial[first(op.bulk_fields)]
  Ω = get_triangulation(Ω)

  m=Gridap.Fields.BlockMap(length(op.bulk_fields),op.bulk_fields)
  uhph_dofs=get_cell_dof_ids(op.trial,Ω)
  # This last step is needed as get_cell_dof_ids(...) returns as many blocks
  # as fields in op.trial, regardless of the FEspaces defined on Ω or not
  uhph_dofs=lazy_map(m,uhph_dofs.args[op.bulk_fields]...)

  uhphₖ=lazy_map(RestrictArrayBlockMap(op.bulk_fields),uhphlhₖ)

  cache = nothing
  free_dof_values=assemble_vector(assem,([lhₑ,uhphₖ],[lhₑ_dofs,uhph_dofs]))
  FEFunction(op.test,free_dof_values), cache
end

# TO-THINK:
# 1. Can TempSkeleton appear more than once? (e.g., the skeleton of different domains)
# 2. Can TempSkeleton appear on the linear form?
# 3. Can BoundaryTriangulation appear on the bilinear form? (e.g., Nitsche BCs)
# 4. May we have a hybridizable weak formulation with triangulations different from TempSkeleton,
#    bulk and BoundaryTriangulation?

# The currently supported scenarios are explicitly encoded in the @check's below.
# These may be modified (along with the code supporting them) as we consider
# more general scenarios.

function _hybridrizable_to_hybrid_contributions(matcontribs,veccontribs)
   Dbi = maximum(map(tr->num_cell_dims(tr), collect(keys(matcontribs.dict))))
   Dli = maximum(map(tr->num_cell_dims(tr), collect(keys(veccontribs.dict))))
   D   = max(Dbi,Dli)

   mskeleton = _find_skeleton(matcontribs); Gridap.Helpers.@check length(mskeleton)==1
   mbulk     = _find_bulk(D,matcontribs)    ; Gridap.Helpers.@check length(mbulk)==1
   mboun     = _find_boundary(matcontribs); Gridap.Helpers.@check length(mboun)==0

   vskeleton = _find_skeleton(veccontribs); Gridap.Helpers.@check length(vskeleton)==0
   vbulk     = _find_bulk(D,veccontribs)    ; Gridap.Helpers.@check length(vbulk)<=1
   vboun     = _find_boundary(veccontribs); Gridap.Helpers.@check length(vboun)<=1

   omatcontribs = DomainContribution()
   oveccontribs = DomainContribution()

   mskeleton = mskeleton[1]
   mbulk     = mbulk[1]

   add_contribution!(omatcontribs,mskeleton,matcontribs[mskeleton])
   add_contribution!(omatcontribs,mskeleton,matcontribs[mbulk])
   if length(vbulk)!=0
    vbulk=vbulk[1]
    add_contribution!(oveccontribs,mskeleton,veccontribs[vbulk])
   end
   if length(vboun)!=0
    add_contribution!(oveccontribs,vboun,veccontribs[vboun])
   end

   omatcontribs, oveccontribs
end

function _find_skeleton(dc::DomainContribution)
  [trian for trian in keys(dc.dict) if isa(trian,TempSkeletonTriangulation)]
end

function _find_bulk(D,dc::DomainContribution)
  [trian for trian in keys(dc.dict)
     if isa(trian,Triangulation{D,D}) && !(isa(trian,TempSkeletonTriangulation))]
end

function _find_boundary(dc::DomainContribution)
  [trian for trian in keys(dc.dict)
     if isa(trian,BoundaryTriangulation)]
end

function _add_static_condensation(matvec,bulk_fields,skeleton_fields)
  Gridap.Helpers.@check length(keys(matvec.dict))==1
  _matvec=DomainContribution()
  for (trian,t) in matvec.dict
    Gridap.Helpers.@check isa(trian,TempSkeletonTriangulation)
    _matvec.dict[trian] = lazy_map(StaticCondensationMap(bulk_fields,skeleton_fields),t)
  end
  _matvec
end

function _block_skeleton_system_contributions(matvec,L::MultiFieldFESpace)
  fdofscw=_get_facet_dofs_cell_wise(L)
  num_cells=length(fdofscw)
  a=fdofscw[1].array
  block_sizes=Fill([length(a[i]) for i=1:length(a)],num_cells)
  m=Scalar2ArrayBlockMap()
  Gridap.Helpers.@check length(keys(matvec.dict))==1
  _matvec=DomainContribution()
  for (trian,t) in matvec.dict
    Gridap.Helpers.@check isa(trian,TempSkeletonTriangulation)
    _matvec.dict[trian] = lazy_map(m,t,block_sizes)
  end
  _matvec
end

function Gridap.FESpaces.get_cell_fe_data(
  fun,
  sface_to_data,
  sglue::FaceToFaceGlue,
  tglue::TempSkeletonGlue)
  model=tglue.trian.model
  cell_wise_facets=_get_cell_wise_facets(model)
  restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets,sface_to_data)
end

function Gridap.FESpaces.get_cell_fe_data(
  fun::typeof(Gridap.FESpaces.get_cell_is_dirichlet),sface_to_data,sglue::FaceToFaceGlue,tglue::TempSkeletonGlue)
  model=tglue.trian.model
  cell_wise_facets=_get_cell_wise_facets(model)
  fdofscb=restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets,sface_to_data)
  _generate_cell_is_dirichlet(fdofscb)
end

function _generate_cell_is_dirichlet(cell_dofs)
  cell_is_dirichlet = Vector{Bool}(undef,length(cell_dofs))
  cache = array_cache(cell_dofs)
  isdirichlet(dof) = dof
  for cell in 1:length(cell_dofs)
    dofs = getindex!(cache,cell_dofs,cell)
    cell_is_dirichlet[cell] = any(isdirichlet,dofs)
  end
  cell_is_dirichlet
end

function _get_facet_dofs_cell_wise(L::Gridap.FESpaces.SingleFieldFESpace)
  Γ=get_triangulation(L)
  model=get_background_model(Γ)
  cell_wise_facets=_get_cell_wise_facets(model)
  restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets,get_cell_dof_ids(L))
end

function _get_facet_dofs_cell_wise(L::MultiFieldFESpace)
  Γ=get_triangulation(L[1])
  model=get_background_model(Γ)
  cell_wise_facets=_get_cell_wise_facets(model)

  facet_fields_dofs_cell_wise=[]
  for S in L
    push!(facet_fields_dofs_cell_wise,
          restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets,get_cell_dof_ids(S)))
  end

  m=BlockMap(num_fields(L),collect(1:num_fields(L)))
  lazy_map(m,facet_fields_dofs_cell_wise...)
end
