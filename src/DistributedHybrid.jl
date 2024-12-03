

# Hybrid distributed model 
function GridapHybrid.Skeleton(model::GridapDistributed.DistributedDiscreteModel)
    cell_ids = get_cell_gids(model)
    trians = map(local_views(model),partition(cell_ids)) do model, indices
        Gridap.Geometry.TriangulationView(
            GridapHybrid.Skeleton(model),
            own_to_local(indices)
        )
    end
    return GridapDistributed.DistributedTriangulation(trians,model)
end

# Integration machinery for distributed hybrid models
const SkeletonView{Dc,Dp} = Gridap.Geometry.TriangulationView{Dc,Dp,<:GridapHybrid.SkeletonTriangulation}

function Gridap.Geometry.is_change_possible(strian::GridapHybrid.SkeletonTriangulation,ttrian::SkeletonView)
    is_change_possible(strian,ttrian.parent)
end

function Gridap.Geometry.is_change_possible(strian::SkeletonView,ttrian::GridapHybrid.SkeletonTriangulation)
    false
end

function Gridap.Geometry.is_change_possible(strian::BodyFittedTriangulation,ttrian::SkeletonView)
    is_change_possible(strian,ttrian.parent)
end

function Gridap.Geometry.is_change_possible(strian::SkeletonView,ttrian::BodyFittedTriangulation)
    false
end

# Keeping ReferenceDomain here due to Julia ambiguity. PhysicalDomain can be incorporated if needed.
function Gridap.CellData.change_domain(
    a::CellField,strian::GridapHybrid.SkeletonTriangulation,::ReferenceDomain,ttrian::SkeletonView,::ReferenceDomain
)
    @assert is_change_possible(strian,ttrian)
    b = change_domain(a,strian,ReferenceDomain(),ttrian.parent,ReferenceDomain())

    data = Gridap.Geometry.restrict(CellData.get_data(b),ttrian.cell_to_parent_cell);
    Gridap.CellData.similar_cell_field(b,data,ttrian,ReferenceDomain())
end

function Gridap.CellData.change_domain(
    a::CellField,strian::BodyFittedTriangulation,::ReferenceDomain,ttrian::SkeletonView,::ReferenceDomain
)
    @assert is_change_possible(strian,ttrian)
    b = change_domain(a,strian,ReferenceDomain(),ttrian.parent,ReferenceDomain())

    data = Gridap.Geometry.restrict(CellData.get_data(b),ttrian.cell_to_parent_cell);
    Gridap.CellData.similar_cell_field(b,data,ttrian,ReferenceDomain())
end

function Gridap.CellData.change_domain(
    a::Gridap.MultiField.MultiFieldFEBasisComponent,
    ttrian::SkeletonView,
    tdomain::DomainStyle)
    if get_triangulation(a) == ttrian # SkeletonView to SkeletonView, return Idem
        return a
    end
    b = change_domain(a, ttrian.parent, tdomain)
    b_sf = b.single_field
    b_data = Gridap.Geometry.restrict(CellData.get_data(b_sf),ttrian.cell_to_parent_cell)
    c_sf = Gridap.CellData.similar_cell_field(b_sf,b_data,ttrian,ReferenceDomain())
    return Gridap.MultiField.MultiFieldFEBasisComponent(b_data,c_sf,a.fieldid,a.nfields)
end

# Cut the data to the local view
function Gridap.FESpaces.get_cell_fe_data(
    fun, f, ttrian::SkeletonView
)
    data = Gridap.FESpaces.get_cell_fe_data(fun,f,ttrian.parent)
    return Gridap.Geometry.restrict(data,ttrian.cell_to_parent_cell)
end

# cellwise normal vector for distributed hybrid models

function GridapHybrid.get_cell_normal_vector(s::GridapDistributed.DistributedTriangulation)
    fields = map(local_views(s)) do si
        parent = get_cell_normal_vector(si.parent)
        data = Geometry.restrict(CellData.get_data(parent),si.cell_to_parent_cell)
        GenericCellField(data,si,ReferenceDomain())
    end
    GridapDistributed.DistributedCellField(fields,s)
end

function GridapHybrid.get_cell_owner_normal_vector(s::GridapDistributed.DistributedTriangulation)
    fields = map(local_views(s)) do si
        parent = get_cell_owner_normal_vector(si.parent)
        data = Geometry.restrict(CellData.get_data(parent),si.cell_to_parent_cell)
        GenericCellField(data,si,ReferenceDomain())
    end
    GridapDistributed.DistributedCellField(fields,s)
end

function GridapHybrid.get_cell_normal_vector(s::GridapHybrid.SkeletonView)
    parent = get_cell_normal_vector(s.parent)
    data = Geometry.restrict(CellData.get_data(parent),s.cell_to_parent_cell)
    GenericCellField(data,s,ReferenceDomain())
end

# For distributed models we need to differentiate between Ωu and Ωp (they have different object ids)
# where u and p are coming from different spaces (example velocity and pressure). The serial implementation 
# does not require this differentiation.
function Geometry.best_target(
    trian1::BodyFittedTriangulation{Dc},
    trian2::BodyFittedTriangulation{Dc}
) where {Dc}
    @check is_change_possible(trian1,trian2)
    @check is_change_possible(trian2,trian1)
    D1 = num_cell_dims(trian1)
    D2 = num_cell_dims(trian2)
    glue1 = get_glue(trian1,Val(D2))
    glue2 = get_glue(trian2,Val(D1))
    best_target(trian1,trian2,glue1,glue2)
end

# Hybrid FE operator for distributed models
function HybridAffineFEOperator(
    weakform::Function,
    trial::GridapDistributed.DistributedMultiFieldFESpace,
    test::GridapDistributed.DistributedMultiFieldFESpace,
    bulk_fields::TB,
    skeleton_fields::TS
) where {TB <: Vector{<:Integer},TS <: Vector{<:Integer}}
    # Invoke weak form of the hybridizable system
    u = get_trial_fe_basis(trial)
    v = get_fe_basis(test)
    biform, liform  = weakform(u, v)
  
    M, L = GridapHybrid._setup_fe_spaces_skeleton_system(trial, test, skeleton_fields)
    
    # Transform DomainContribution objects of the hybridizable system into a
    # suitable form for assembling the linear system defined on the skeleton
    # (i.e., the hybrid system)
    data = map(local_views(biform), local_views(liform),local_views(M),local_views(L)) do biform, liform, M, L
        obiform, oliform = GridapHybrid._merge_bulk_and_skeleton_contributions(biform, liform)
        # Pair LHS and RHS terms associated to SkeletonTriangulation
        matvec, mat, vec = Gridap.FESpaces._pair_contribution_when_possible(obiform, oliform)
        # Add StaticCondensationMap to matvec terms
        matvec = GridapHybrid._add_static_condensation(matvec, bulk_fields, skeleton_fields)
  
        if (length(skeleton_fields) != 1)
            matvec = GridapHybrid._block_skeleton_system_contributions(matvec, L)
        end

        uhd = zero(M)
        matvec, mat = Gridap.FESpaces._attach_dirichlet(matvec, mat, uhd)

        data = Gridap.FESpaces._collect_cell_matrix_and_vector(M, L, matvec, mat, vec)
        return data
    end 
  
    assem = SparseMatrixAssembler(M, L)
    A, b = assemble_matrix_and_vector(assem, data)

    skeleton_op = AffineFEOperator(M, L, A, b)
    HybridAffineFEOperator(weakform, trial, test, bulk_fields, skeleton_fields, skeleton_op)
end
  

# solve function for hybrid distributed models
function Gridap.FESpaces.solve!(uh::GridapDistributed.DistributedMultiFieldCellField, solver::LinearFESolver, op::HybridAffineFEOperator, cache)
    # Solve linear system defined on the skeleton
    lh = solve(op.skeleton_op)
  
    # Invoke weak form of the hybridizable system
    u = get_trial_fe_basis(op.trial)
    v = get_fe_basis(op.test)
    biform, liform  = op.weakform(u, v)
  
    # Transform DomainContribution objects of the hybridizable system into a
    # suitable form for assembling the linear system defined on the skeleton
    # (i.e., the hybrid system) and pair LHS and RHS terms associated to 
    # SkeletonTriangulation
    matvec = map(local_views(biform), local_views(liform)) do biformi, liformi
      obiform, oliform = GridapHybrid._merge_bulk_and_skeleton_contributions(biformi, liformi)
      matvec, _, _ = Gridap.FESpaces._pair_contribution_when_possible(obiform, oliform)
      return matvec
    end
  
    free_dof_values = GridapHybrid._compute_hybridizable_from_skeleton_free_dof_values(
        lh,
        op.trial,
        op.test,
        local_views(Gridap.FESpaces.get_trial(op.skeleton_op)),
        matvec,
        op.bulk_fields,
        op.skeleton_fields
    )
      
    ids = partition(get_free_dof_ids(op.trial))
    xh  =  FEFunction(op.trial, PVector(free_dof_values,ids))              
    cache = nothing
    return xh, cache
end

function GridapHybrid._compute_hybridizable_from_skeleton_free_dof_values(
    skeleton_fe_function::GridapDistributed.DistributedCellField,
    trial_hybridizable::GridapDistributed.DistributedMultiFieldFESpace,
    test_hybridizable::GridapDistributed.DistributedMultiFieldFESpace,
    trial_skeleton,
    matvec::AbstractArray{<:DomainContribution},
    bulk_fields::Vector{Int},
    skeleton_fields::Vector{Int}
)
    free_dof_values = map(local_views(skeleton_fe_function), local_views(trial_hybridizable), local_views(test_hybridizable), local_views(matvec)) do skeleton_fe_function, trial_hybridizable, test_hybridizable, matvec
        # Convert dof-wise dof values of lh into cell-wise dof values lhₖ
        Γ = first(keys(matvec.dict))
        Gridap.Helpers.@check isa(Γ, GridapHybrid.SkeletonTriangulation) || isa(Γ, GridapHybrid.SkeletonView)
        lhₖ = get_cell_dof_values(skeleton_fe_function, Γ)
        
        # Compute cell-wise dof values of bulk fields out of lhₖ
        t = matvec.dict[Γ]
        m = BackwardStaticCondensationMap(bulk_fields, skeleton_fields)
        uhphlhₖ = lazy_map(m, t, lhₖ)
        
        model = get_background_model(Γ)
        cell_wise_facets_parent = GridapHybrid._get_cell_wise_facets(model)
        cell_wise_facets = Gridap.Geometry.restrict(cell_wise_facets_parent, Γ.cell_to_parent_cell)
        cells_around_facets_parent = GridapHybrid._get_cells_around_facets(model)
        facets = unique(vcat(cell_wise_facets...))
        
        nfields = length(bulk_fields) + length(skeleton_fields)
        m = Gridap.Fields.BlockMap(nfields, skeleton_fields)
        L = Gridap.FESpaces.get_fe_space(skeleton_fe_function)
        lhₑ = lazy_map(m,
                        GridapHybrid.convert_cell_wise_dofs_array_to_facet_dofs_array(
                        cells_around_facets_parent,
                        cell_wise_facets_parent,
                        lhₖ,
                        get_cell_dof_ids(L))...)
        lhₑ = Gridap.Geometry.restrict(lhₑ, 1:length(facets))            
        
        lhₑ_dofs = Gridap.FESpaces.get_cell_dof_ids(trial_hybridizable, get_triangulation(L))
        lhₑ_dofs = lazy_map(m, lhₑ_dofs.args[skeleton_fields]...)
        lhₑ_dofs = Gridap.Geometry.restrict(lhₑ_dofs, facets)
        
        Ω = trial_hybridizable[first(bulk_fields)]
        Ω = get_triangulation(Ω)
        
        m = Gridap.Fields.BlockMap(length(bulk_fields), bulk_fields)
        uhph_dofs = get_cell_dof_ids(trial_hybridizable, Ω)
        
        # This last step is needed as get_cell_dof_ids(...) returns as many blocks
        # as fields in trial, regardless of the FEspaces defined on Ω or not
        uhph_dofs = lazy_map(m, uhph_dofs.args[bulk_fields]...)
        uhph_dofs = Gridap.Geometry.restrict(uhph_dofs,Γ.cell_to_parent_cell)
        
        uhphₖ = lazy_map(RestrictArrayBlockMap(bulk_fields), uhphlhₖ)             
        vecdata = ([lhₑ,uhphₖ], [lhₑ_dofs,uhph_dofs])
                        
        assem = SparseMatrixAssembler(trial_hybridizable, test_hybridizable)
        free_dof_values = assemble_vector(assem, vecdata)       
        return free_dof_values          
    end    
end

# HHO specific functionalities

function (op::AbstractArray{<:ReconstructionFEOperator})(v::GridapDistributed.DistributedMultiFieldCellField) 
    image_span = map(local_views(op), local_views(v)) do opi, vi
      basis_style = GridapHybrid._get_basis_style(vi)
      LHSf, RHSf = GridapHybrid._evaluate_forms(opi, vi)
      cell_dofs = GridapHybrid._compute_cell_dofs(LHSf, RHSf)
      O = opi.test_space
      GridapHybrid._generate_image_space_span(opi,O,vi,cell_dofs,basis_style)
    end
    return image_span
end

function (op::AbstractArray{<:ProjectionFEOperator})(v::GridapDistributed.DistributedMultiFieldCellField)
    image_span = map(local_views(op), local_views(v)) do opi, vi
        basis_style = GridapHybrid._get_basis_style(vi)
        LHSf,RHSf= GridapHybrid._evaluate_forms(opi,vi)
        cell_dofs= GridapHybrid._compute_cell_dofs(LHSf,RHSf)
        # We get the test space so that we don't deal with
        # the transpose underlying the trial basis
        O = opi.test_space
        GridapHybrid._generate_image_space_span(opi,O,vi,cell_dofs,basis_style)
    end
    return image_span
end
