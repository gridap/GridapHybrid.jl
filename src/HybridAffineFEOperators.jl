struct HybridAffineFEOperator{TB,TS} <: FEOperator
    weakform::Function
    trial::MultiFieldFESpace
    test::MultiFieldFESpace
    bulk_fields::TB
    skeleton_fields::TS
    skeleton_op::AffineFEOperator
end

function HybridAffineFEOperator(
  weakform::Function,
  trial::MultiFieldFESpace,
  test::MultiFieldFESpace,
  bulk_fields::TB,
  skeleton_fields::TS) where {TB <: Vector{<:Integer},TS <: Vector{<:Integer}}

  # Invoke weak form of the hybridizable system
    u = get_trial_fe_basis(trial)
    v = get_fe_basis(test)
    biform, liform  = weakform(u, v)

  # Transform DomainContribution objects of the hybridizable system into a
  # suitable form for assembling the linear system defined on the skeleton
  # (i.e., the hybrid system)
    obiform, oliform = _merge_bulk_and_skeleton_contributions(biform, liform)

  # Pair LHS and RHS terms associated to SkeletonTriangulation
    matvec, mat, vec = Gridap.FESpaces._pair_contribution_when_possible(obiform, oliform)

  # Add StaticCondensationMap to matvec terms
    matvec = _add_static_condensation(matvec, bulk_fields, skeleton_fields)

    M, L = _setup_fe_spaces_skeleton_system(trial, test, skeleton_fields)

    if (length(skeleton_fields) != 1)
        matvec = _block_skeleton_system_contributions(matvec, L)
    end
    assem = SparseMatrixAssembler(M, L)

  # Attach strong imposition of Dirichlet BCs
    uhd = zero(M)
    matvec, mat = Gridap.FESpaces._attach_dirichlet(matvec, mat, uhd)

    data = Gridap.FESpaces._collect_cell_matrix_and_vector(M, L, matvec, mat, vec)

    A, b = assemble_matrix_and_vector(assem, data)

    skeleton_op = AffineFEOperator(M, L, A, b)
    HybridAffineFEOperator(weakform, trial, test, bulk_fields, skeleton_fields, skeleton_op)
end

function _setup_fe_spaces_skeleton_system(trial, test, skeleton_fields)
    if (length(skeleton_fields) == 1)
    # Single-field skeleton system
        M = trial[skeleton_fields[1]]
        L = test[skeleton_fields[1]]
    else
    # Multi-field skeleton system
        M = MultiFieldFESpace([trial[i] for i in skeleton_fields])
        L = MultiFieldFESpace([test[i] for i in skeleton_fields])
    end
    M, L
end

Gridap.FESpaces.get_test(feop::HybridAffineFEOperator) = feop.test
Gridap.FESpaces.get_trial(feop::HybridAffineFEOperator) = feop.trial
function Gridap.FESpaces.solve(op::HybridAffineFEOperator)
    solver = LinearFESolver()
    solve(solver, op)
end

function Gridap.FESpaces.solve!(uh, solver::LinearFESolver, op::HybridAffineFEOperator, cache)
   # Solve linear system defined on the skeleton.
   # Some solvers do not support 0x0 matrices. Avoid
   # calling solve in such an scenario.
   if (size(op.skeleton_op.op.matrix) != (0,0))
    lh = solve(op.skeleton_op)
   else
    lh = FEFunction(op.skeleton_op.trial,zeros(num_free_dofs(op.skeleton_op.trial)))
   end

  # Invoke weak form of the hybridizable system
    u = get_trial_fe_basis(op.trial)
    v = get_fe_basis(op.test)
    biform, liform  = op.weakform(u, v)

  # Transform DomainContribution objects of the hybridizable system into a
  # suitable form for assembling the linear system defined on the skeleton
  # (i.e., the hybrid system)
    obiform, oliform = _merge_bulk_and_skeleton_contributions(biform, liform)

  # Pair LHS and RHS terms associated to SkeletonTriangulation
    matvec, _, _ = Gridap.FESpaces._pair_contribution_when_possible(obiform, oliform)

    free_dof_values = _compute_hybridizable_from_skeleton_free_dof_values(
                    lh,
                    op.trial,
                    op.test,
                    Gridap.FESpaces.get_trial(op.skeleton_op),
                    matvec,
                    op.bulk_fields,
                    op.skeleton_fields)

    cache = nothing
    FEFunction(op.trial, free_dof_values), cache
end

function _compute_hybridizable_from_skeleton_free_dof_values(skeleton_fe_function,
                                                             trial_hybridizable,
                                                             test_hybridizable,
                                                             trial_skeleton,
                                                             matvec,
                                                             bulk_fields,
                                                             skeleton_fields)

  # Convert dof-wise dof values of lh into cell-wise dof values lhₖ
    Γ = first(keys(matvec.dict))
    Gridap.Helpers.@check isa(Γ, SkeletonTriangulation)
    lhₖ = get_cell_dof_values(skeleton_fe_function, Γ)

  # Compute cell-wise dof values of bulk fields out of lhₖ
    t = matvec.dict[Γ]
    m = BackwardStaticCondensationMap(bulk_fields, skeleton_fields)
    uhphlhₖ = lazy_map(m, t, lhₖ)

    model = get_background_model(Γ)
    cell_wise_facets = _get_cell_wise_facets(model)
    cells_around_facets = _get_cells_around_facets(model)

    nfields = length(bulk_fields) + length(skeleton_fields)
    m = Gridap.Fields.BlockMap(nfields, skeleton_fields)
    L = Gridap.FESpaces.get_fe_space(skeleton_fe_function)
    lhₑ = lazy_map(m,
                 convert_cell_wise_dofs_array_to_facet_dofs_array(
                 cells_around_facets,
                 cell_wise_facets,
                 lhₖ,
                 get_cell_dof_ids(L))...)

    assem = SparseMatrixAssembler(trial_hybridizable, test_hybridizable)
    lhₑ_dofs = get_cell_dof_ids(trial_hybridizable, get_triangulation(L))
    lhₑ_dofs = lazy_map(m, lhₑ_dofs.args[skeleton_fields]...)

    Ω = trial_hybridizable[first(bulk_fields)]
    Ω = get_triangulation(Ω)

    m = Gridap.Fields.BlockMap(length(bulk_fields), bulk_fields)
    uhph_dofs = get_cell_dof_ids(trial_hybridizable, Ω)
  # This last step is needed as get_cell_dof_ids(...) returns as many blocks
  # as fields in trial, regardless of the FEspaces defined on Ω or not
    uhph_dofs = lazy_map(m, uhph_dofs.args[bulk_fields]...)

    uhphₖ = lazy_map(RestrictArrayBlockMap(bulk_fields), uhphlhₖ)

    free_dof_values = assemble_vector(assem, ([lhₑ,uhphₖ], [lhₑ_dofs,uhph_dofs]))
end

function _get_cell_wise_facets(model::DiscreteModel)
    gtopo = get_grid_topology(model)
    D     = num_cell_dims(model)
    Gridap.Geometry.get_faces(gtopo, D, D - 1)
end

function _get_cells_around_facets(model::DiscreteModel)
    gtopo = get_grid_topology(model)
    D     = num_cell_dims(model)
    Gridap.Geometry.get_faces(gtopo, D - 1, D)
end

struct ExtractFacetDofsFromCellDofsMap{T <: AbstractVector{<:AbstractArray}} <: Gridap.Fields.Map
    cell_dofs::T
end

function Gridap.Arrays.return_cache(
       k::ExtractFacetDofsFromCellDofsMap{<:AbstractVector{<:AbstractVector}},
       cellgid_facetlpos_ndofs::NTuple{3})
    cell_dofs_cache  = Gridap.Arrays.array_cache(k.cell_dofs)
    T = eltype(eltype(k.cell_dofs))
    facet_dofs_cache = Gridap.Arrays.CachedArray(zeros(T, cellgid_facetlpos_ndofs[3]))
    cell_dofs_cache, facet_dofs_cache
end

function Gridap.Arrays.evaluate!(cache,
       k::ExtractFacetDofsFromCellDofsMap{<:AbstractVector{<:AbstractVector}},
       cellgid_facetlpos_ndofs::NTuple{3})
    cell_dofs_cache, facet_dofs_cache = cache
    cellgid, facetlpos, ndofs = cellgid_facetlpos_ndofs
    Gridap.Arrays.setsize!(facet_dofs_cache, (ndofs,))
    facet_dofs  = facet_dofs_cache.array
    cell_dofs   = Gridap.Arrays.getindex!(cell_dofs_cache, k.cell_dofs, cellgid)
    facet_dofs .= cell_dofs[facetlpos:facetlpos + ndofs - 1]
end

function Gridap.Arrays.return_cache(
       k::ExtractFacetDofsFromCellDofsMap{<:AbstractVector{<:AbstractMatrix}},
       cellgid_facetlpos_ndofs::NTuple{3})
    cell_dofs_cache  = Gridap.Arrays.array_cache(k.cell_dofs)
    cell_dofs        = getindex!(cell_dofs_cache,k.cell_dofs,1)
    T = eltype(eltype(k.cell_dofs))
    facet_dofs_cache = Gridap.Arrays.CachedArray(zeros(T,
                                                 cellgid_facetlpos_ndofs[3], size(cell_dofs,2)))
    cell_dofs_cache, facet_dofs_cache
end

function Gridap.Arrays.evaluate!(cache,
       k::ExtractFacetDofsFromCellDofsMap{<:AbstractVector{<:AbstractMatrix}},
       cellgid_facetlpos_ndofs::NTuple{3})
    cell_dofs_cache, facet_dofs_cache = cache
    cellgid, facetlpos, ndofs = cellgid_facetlpos_ndofs
    cell_dofs   = Gridap.Arrays.getindex!(cell_dofs_cache, k.cell_dofs, cellgid)
    Gridap.Arrays.setsize!(facet_dofs_cache, (ndofs,size(cell_dofs,2)))
    facet_dofs  = facet_dofs_cache.array
    facet_dofs .= cell_dofs[facetlpos:facetlpos + ndofs - 1,:]
end

function convert_cell_wise_dofs_array_to_facet_dofs_array(
      cells_around_facets,
      cell_wise_facets,
      cell_dofs::AbstractVector{<:AbstractArray{<:Real}},
      facet_dofs_ids::AbstractVector{<:AbstractVector{<:Integer}})
    glue = _generate_glue_among_facet_and_cell_wise_dofs_arrays(
   cells_around_facets, cell_wise_facets, facet_dofs_ids)
    k = ExtractFacetDofsFromCellDofsMap(cell_dofs)
    [lazy_map(k, glue)]
end

function _generate_glue_among_facet_and_cell_wise_dofs_arrays(
  cells_around_facets,
  cell_wise_facets,
  facet_dof_ids::AbstractVector{<:AbstractVector{<:Integer}})

    c1 = array_cache(cells_around_facets)
    c2 = array_cache(cell_wise_facets)
    c3 = array_cache(facet_dof_ids)

    result = Vector{NTuple{3,Int}}(undef, length(facet_dof_ids))
    current = 1
    ndofs = 0
    for facet_gid = 1:length(cells_around_facets)
        cell_gid = Gridap.Arrays.getindex!(c1, cells_around_facets, facet_gid)[1]
        current_cell_facets = Gridap.Arrays.getindex!(c2, cell_wise_facets, cell_gid)
        pos = 1
        for facet_gid_in_cell in current_cell_facets
            ndofs = length(Gridap.Arrays.getindex!(c3, facet_dof_ids, facet_gid_in_cell))
            if (facet_gid == facet_gid_in_cell)
                break
            else
                pos = pos + ndofs
            end
        end
        result[facet_gid] = (cell_gid, pos, ndofs)
    end
    result
end

function convert_cell_wise_dofs_array_to_facet_dofs_array(
 cells_around_facets,
 cell_wise_facets,
 cell_dofs::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}},
 facet_dofs_ids::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}})
    Gridap.Helpers.@check cell_dofs.maps.value.size == facet_dofs_ids.maps.value.size

    [convert_cell_wise_dofs_array_to_facet_dofs_array(
     cells_around_facets,
     cell_wise_facets,
     field_cell_dofs,
     field_facet_dof_ids)[1] for (field_cell_dofs, field_facet_dof_ids) in
                                  zip(cell_dofs.args, facet_dofs_ids.args)]
end


# TO-THINK:
# 1. Can Skeleton appear more than once? (e.g., the skeleton of different domains)
# 2. Can Skeleton appear on the linear form?
# 3. Can BoundaryTriangulation appear on the bilinear form? (e.g., Nitsche BCs)
# 4. May we have a hybridizable weak formulation with triangulations different from Skeleton,
#    bulk and BoundaryTriangulation?

# The currently supported scenarios are explicitly encoded in the @check's below.
# These may be modified (along with the code supporting them) as we consider
# more general scenarios.

function _merge_bulk_and_skeleton_contributions(matcontribs, veccontribs)
    Dbi = maximum(map(tr -> num_cell_dims(tr), collect(keys(matcontribs.dict))))
    Dli = maximum(map(tr -> num_cell_dims(tr), collect(keys(veccontribs.dict))))
    D   = max(Dbi, Dli)

    mskeleton = _find_skeleton(matcontribs); Gridap.Helpers.@check length(mskeleton) == 1
    mbulk     = _find_bulk(D, matcontribs)    ; Gridap.Helpers.@check length(mbulk) <= 1
    mboun     = _find_boundary(matcontribs); Gridap.Helpers.@check length(mboun) == 0

    vskeleton = _find_skeleton(veccontribs); Gridap.Helpers.@check length(vskeleton) <= 1
    vbulk     = _find_bulk(D, veccontribs)    ; Gridap.Helpers.@check length(vbulk) <= 1
    vboun     = _find_boundary(veccontribs); Gridap.Helpers.@check length(vboun) <= 1

    omatcontribs = DomainContribution()
    oveccontribs = DomainContribution()

    if length(mbulk) != 0
        add_contribution!(omatcontribs, mskeleton..., matcontribs[mbulk...])
    end
    add_contribution!(omatcontribs, mskeleton..., matcontribs[mskeleton...])

    if length(vbulk) != 0
        add_contribution!(oveccontribs, mskeleton..., veccontribs[vbulk...])
    end
    if length(vskeleton) != 0
        add_contribution!(oveccontribs, vskeleton..., veccontribs[vskeleton...])
    end
    if length(vboun) != 0
        add_contribution!(oveccontribs, vboun..., veccontribs[vboun...])
    end

    omatcontribs, oveccontribs
end

function _merge_bulk_and_skeleton_vector_contributions(veccontribs)
    Dli = maximum(map(tr -> num_cell_dims(tr), collect(keys(veccontribs.dict))))
    D   = Dli
    vskeleton = _find_skeleton(veccontribs); Gridap.Helpers.@check length(vskeleton) == 1
    vbulk     = _find_bulk(D, veccontribs)    ; Gridap.Helpers.@check length(vbulk) <= 1
    vboun     = _find_boundary(veccontribs); Gridap.Helpers.@check length(vboun) <= 1
    oveccontribs = DomainContribution()
    if length(vbulk) != 0
        vbulk = vbulk[1]
        add_contribution!(oveccontribs, vskeleton..., veccontribs[vbulk])
    end
    if length(vboun) != 0
        add_contribution!(oveccontribs, vboun..., veccontribs[vboun])
    end
    oveccontribs
end

function _merge_bulk_and_skeleton_matrix_contributions(matcontribs)
    Dbi = maximum(map(tr -> num_cell_dims(tr), collect(keys(matcontribs.dict))))
    D   = Dbi
    mskeleton = _find_skeleton(matcontribs); Gridap.Helpers.@check length(mskeleton) == 1
    mbulk     = _find_bulk(D, matcontribs)    ; Gridap.Helpers.@check length(mbulk) == 1
    mboun     = _find_boundary(matcontribs); Gridap.Helpers.@check length(mboun) == 0
    omatcontribs = DomainContribution()
    mskeleton = mskeleton[1]
    mbulk     = mbulk[1]
    add_contribution!(omatcontribs, mskeleton, matcontribs[mskeleton])
    add_contribution!(omatcontribs, mskeleton, matcontribs[mbulk])
    omatcontribs
end

function _find_skeleton(dc::DomainContribution)
    [trian for trian in keys(dc.dict) if isa(trian, SkeletonTriangulation)]
end

function _find_bulk(D, dc::DomainContribution)
    [trian for trian in keys(dc.dict)
     if isa(trian, Triangulation{D}) && !(isa(trian, SkeletonTriangulation))]
end

function _find_boundary(dc::DomainContribution)
    [trian for trian in keys(dc.dict)
     if isa(trian, BoundaryTriangulation)]
end

function _add_static_condensation(matvec, bulk_fields, skeleton_fields)
    Gridap.Helpers.@check length(keys(matvec.dict)) == 1
    _matvec = DomainContribution()
    for (trian, t) in matvec.dict
        Gridap.Helpers.@check isa(trian, SkeletonTriangulation)
        _matvec.dict[trian] = lazy_map(StaticCondensationMap(bulk_fields, skeleton_fields), t)
    end
    _matvec
end

function _block_skeleton_system_contributions(matvec, L::MultiFieldFESpace)
    fdofscw = _get_facet_dofs_cell_wise(L)
    num_cells = length(fdofscw)
    a = fdofscw[1].array
    block_sizes = Fill([length(a[i]) for i = 1:length(a)], num_cells)
    m = Scalar2ArrayBlockMap()
    Gridap.Helpers.@check length(keys(matvec.dict)) == 1
    _matvec = DomainContribution()
    for (trian, t) in matvec.dict
        Gridap.Helpers.@check isa(trian, SkeletonTriangulation)
        _matvec.dict[trian] = lazy_map(m, t, block_sizes)
    end
    _matvec
end

function Gridap.FESpaces.get_cell_fe_data(
  fun,
  sface_to_data,
  sglue::FaceToFaceGlue{Dc},
  tglue::SkeletonGlue) where Dc
    model = tglue.trian.model
    if Dc == num_cell_dims(model)
        sface_to_data
    else
        @assert Dc == num_cell_dims(model) - 1
        cell_wise_facets = _get_cell_wise_facets(model)
        restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets, sface_to_data)
    end
end

# TO-DO: We need additional expressivity in Gridap
# This is required in order to be able to evaluate the residual
# of the hybridizable formulation. Commented out while we do not
# perform the necessary modifications to Gridap

# D-dim strian
# No need to transform sface_to_data
# function Gridap.FESpaces.get_cell_fe_data(
#   fun::typeof(get_cell_dof_ids),
#   sface_to_data,
#   sglue::FaceToFaceGlue{<:Gridap.Arrays.IdentityVector},
#   tglue::SkeletonGlue)
#   sface_to_data
# end

struct RestrictFacetDoFsToSkeleton{F} <: Gridap.Fields.Map
    facet_dofs::F
end

function Gridap.Arrays.return_cache(
  k::RestrictFacetDoFsToSkeleton,
  cell_facets::AbstractArray{<:Integer})
    if (isa(k.facet_dofs, AbstractVector{<:Number}))
        T = eltype(k.facet_dofs)
    elseif (isa(k.facet_dofs, AbstractVector{<:AbstractVector{<:Number}}))
    T = eltype(eltype(k.facet_dofs))
  else
    Gridap.Helpers.@check false
    end
    array_cache(k.facet_dofs), Gridap.Arrays.CachedVector(T)
end

function Gridap.Arrays.evaluate!(
  cache,
  k::RestrictFacetDoFsToSkeleton,
  cell_facets::AbstractArray{<:Integer})
    fdc, c = cache
    Gridap.Arrays.setsize!(c, (_count_dofs_cell(cell_facets, k.facet_dofs, fdc),))
    _fill_dofs_cell!(c.array, cell_facets, k.facet_dofs, fdc)
    c.array
end

function _count_dofs_cell(cell_facets, facet_dofs, facet_dofs_cache)
    count = 0
    for facet_id in cell_facets
        current_facet_dofs = getindex!(facet_dofs_cache, facet_dofs, facet_id)
        count = count + length(current_facet_dofs)
    end
    count
end

function _fill_dofs_cell!(cell_dofs, cell_facets, facet_dofs, facet_dofs_cache)
    current = 1
    for facet_id in cell_facets
        current_facet_dofs = getindex!(facet_dofs_cache, facet_dofs, facet_id)
        for i in current_facet_dofs
            cell_dofs[current] = i
            current = current + 1
        end
    end
    current
end

function restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets, facet_dof_ids)
    m = RestrictFacetDoFsToSkeleton(facet_dof_ids)
    lazy_map(m, cell_wise_facets)
end

function Gridap.FESpaces.get_cell_fe_data(
  fun::typeof(Gridap.FESpaces.get_cell_is_dirichlet),sface_to_data,sglue::FaceToFaceGlue,tglue::SkeletonGlue)
    model = tglue.trian.model
    cell_wise_facets = _get_cell_wise_facets(model)
    fdofscb = restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets, sface_to_data)
    _generate_cell_is_dirichlet(fdofscb)
end

function _generate_cell_is_dirichlet(cell_dofs)
    cell_is_dirichlet = Vector{Bool}(undef, length(cell_dofs))
    cache = array_cache(cell_dofs)
    isdirichlet(dof) = dof
    for cell in 1:length(cell_dofs)
        dofs = getindex!(cache, cell_dofs, cell)
        cell_is_dirichlet[cell] = any(isdirichlet, dofs)
    end
    cell_is_dirichlet
end

function _get_facet_dofs_cell_wise(L::Gridap.FESpaces.SingleFieldFESpace)
    Γ = get_triangulation(L)
    model = get_background_model(Γ)
    cell_wise_facets = _get_cell_wise_facets(model)
    restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets, get_cell_dof_ids(L))
end

function _get_facet_dofs_cell_wise(L::MultiFieldFESpace)
    Γ = get_triangulation(L[1])
    model = get_background_model(Γ)
    cell_wise_facets = _get_cell_wise_facets(model)

    facet_fields_dofs_cell_wise = []
    for S in L
        push!(facet_fields_dofs_cell_wise,
          restrict_facet_dof_ids_to_cell_boundary(cell_wise_facets, get_cell_dof_ids(S)))
    end

    m = BlockMap(num_fields(L), collect(1:num_fields(L)))
    lazy_map(m, facet_fields_dofs_cell_wise...)
end
