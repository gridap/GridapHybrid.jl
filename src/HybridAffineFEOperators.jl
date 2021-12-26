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

  # TO-DO: get rid of this temporary limitation
  @assert length(skeleton_fields)==1

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

  M = trial[skeleton_fields[1]]
  L = test[skeleton_fields[1]]
  #@assert length(M)==1
  #@assert length(L)==1
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

  # Convert dof-wise dof values of lh into cell-wise dof values lhₖ
  L=Gridap.FESpaces.get_fe_space(lh)
  Γ=get_triangulation(L)
  model=get_background_model(Γ)
  cell_wise_facets=_get_cell_wise_facets(model)
  fdofscb=restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets,get_cell_dof_ids(L))
  m=Broadcasting(Gridap.Fields.PosNegReindex(
                 Gridap.FESpaces.get_free_dof_values(lh),
                 lh.dirichlet_values))
  lhₖ= lazy_map(m,fdofscb)

  # Compute cell-wise dof values of bulk fields out of lhₖ

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

  Γ=first(keys(matvec.dict))
  @assert isa(Γ,TempSkeletonTriangulation)
  t = matvec.dict[Γ]
  m=BackwardStaticCondensationMap(op.bulk_fields,op.skeleton_fields)
  uhphlhₖ=lazy_map(m,t,lhₖ)

  cell_wise_facets=_get_cell_wise_facets(model)
  cells_around_facets=_get_cells_around_facets(model)

  nfields=length(op.bulk_fields)+length(op.skeleton_fields)
  ### To consider first(op.skeleton_fields)
  lhₑ=lazy_map(Gridap.Fields.BlockMap(nfields,first(op.skeleton_fields)),ExploringGridapHybridization.convert_cell_wise_dofs_array_to_facet_dofs_array(
      cells_around_facets,
      cell_wise_facets,
      lhₖ,
      get_cell_dof_ids(L)))

  assem = SparseMatrixAssembler(op.trial,op.test)
  lhₑ_dofs=get_cell_dof_ids(op.trial,get_triangulation(L))

  Ω = op.trial[first(op.bulk_fields)]
  Ω = get_triangulation(Ω)

  m=Gridap.Fields.BlockMap(length(op.bulk_fields),op.bulk_fields)
  uhph_dofs=get_cell_dof_ids(op.trial,Ω)
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

# The currently supported scenarios are explicitly encoded in the @assert's below.
# These may be modified (along with the code supporting them) as we consider
# more general scenarios.

function _hybridrizable_to_hybrid_contributions(matcontribs,veccontribs)
   Dbi = maximum(map(tr->num_cell_dims(tr), collect(keys(matcontribs.dict))))
   Dli = maximum(map(tr->num_cell_dims(tr), collect(keys(veccontribs.dict))))
   D   = max(Dbi,Dli)

   mskeleton = _find_skeleton(matcontribs); @assert length(mskeleton)==1
   mbulk     = _find_bulk(D,matcontribs)    ; @assert length(mbulk)==1
   mboun     = _find_boundary(matcontribs); @assert length(mboun)==0

   vskeleton = _find_skeleton(veccontribs); @assert length(vskeleton)==0
   vbulk     = _find_bulk(D,veccontribs)    ; @assert length(vbulk)<=1
   vboun     = _find_boundary(veccontribs); @assert length(vboun)<=1

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
  @assert length(keys(matvec.dict))==1
  _matvec=DomainContribution()
  for (trian,t) in matvec.dict
    @assert isa(trian,TempSkeletonTriangulation)
    _matvec.dict[trian] = lazy_map(StaticCondensationMap(bulk_fields,skeleton_fields),t)
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
