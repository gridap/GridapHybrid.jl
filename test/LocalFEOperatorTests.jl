module LocalFEOperatorTests
  using Gridap
  using GridapHybrid
  using Test
  using LinearAlgebra

  function setup_reduction_operator(UK_U∂K,VK_V∂K,dΩ,d∂K)
    m( (uK,u∂K), (vK,v∂K) ) = ∫(vK*uK)dΩ  + ∫(v∂K*u∂K)d∂K
    n( uhK     , (vK,v∂K) ) = ∫(vK*uhK)dΩ + ∫(v∂K*uhK)d∂K
    LocalFEOperator((m,n),UK_U∂K,VK_V∂K)
  end

  function setup_l2_projection_operator(UK_U∂K,VK_V∂K,dΩ,d∂K)
    m((uK,u∂K)  , (vK,v∂K)) = ∫(vK*uK)dΩ  + ∫(v∂K*u∂K)d∂K
    n((uhK,uh∂K), (vK,v∂K)) = ∫(vK*uhK)dΩ + ∫(v∂K*uh∂K)d∂K
    LocalFEOperator((m,n),UK_U∂K,VK_V∂K)
  end

  function setup_difference_operator(l2_projection_op,reconstruction_op)
     function _op(uK_u∂K)
        u_rec=reconstruction_op(uK_u∂K)
        uhK,uh∂K=uK_u∂K
        l2_projection_op((u_rec-uhK,u_rec-uh∂K))
     end
     return _op
  end

  function setup_skeleton_difference_operator(reduction_op,reconstruction_op)
    function _op(uK_u∂K)
       uK,u∂K=uK_u∂K
       u_rec=reconstruction_op(uK_u∂K)
       reduction_op((u_rec-uK))
    end
    return _op
 end

  function setup_reconstruction_operator(model, order, dΩ, d∂K)
    T         = Float64
    nK        = get_cell_normal_vector(d∂K.quad.trian)
    refferecᵤ = ReferenceFE(orthogonal_basis, T, order+1)
    reffe_nzm = ReferenceFE(orthogonal_basis, T, order+1; subspace=:NonZeroMean)
    reffe_zm  = ReferenceFE(orthogonal_basis, T, order+1; subspace=:ZeroMean)
    reffe_c   = ReferenceFE(monomial_basis  , T, order+1; subspace=:OnlyConstant)
    reffe_nc  = ReferenceFE(monomial_basis  , T, order+1; subspace=:ExcludeConstant)

    VKR     = TestFESpace(Ω  , refferecᵤ; conformity=:L2)
    UKR     = TrialFESpace(VKR)

    UKR_NZM = TrialFESpace(TestFESpace(Ω, reffe_nzm; conformity=:L2))
    UKR_ZM  = TrialFESpace(TestFESpace(Ω, reffe_zm; conformity=:L2))
    VKR_C   = TestFESpace(Ω, reffe_c ; conformity=:L2, vector_type=Vector{Float64})
    VKR_NC  = TestFESpace(Ω, reffe_nc; conformity=:L2, vector_type=Vector{Float64})

    VKR_DS_DECOMP = MultiFieldFESpace([VKR_C,VKR_NC])
    UKR_DS_DECOMP = MultiFieldFESpace([UKR_NZM,UKR_ZM])


    m( (u_nzm,u_zm), (v_c,v_nc) ) = ∫(∇(v_nc)⋅∇(u_zm))dΩ +
                                    ∫(∇(v_nc)⋅∇(u_nzm))dΩ +
                                    ∫(v_c*u_nzm)dΩ
    n( (uK,u∂K), (v_c,v_nc)     ) = ∫(-Δ(v_nc)*uK)dΩ + ∫(v_c*uK)dΩ + ∫((∇(v_nc)⋅nK)*u∂K)d∂K

    LocalFEOperator((m,n),UKR,VKR;
                    trial_space_ds_decomp=UKR_DS_DECOMP,
                    test_space_ds_decomp=VKR_DS_DECOMP)
  end

  model=CartesianDiscreteModel((0,1,0,1),(2,2))

  # Geometry
  D  = num_cell_dims(model)
  Ω  = Triangulation(ReferenceFE{D},model)
  Γ  = Triangulation(ReferenceFE{D-1},model)
  ∂K = GridapHybrid.Skeleton(model)

  # Reference FEs
  order     = 1
  reffeᵤ    = ReferenceFE(lagrangian,Float64,order  ;space=:P)

  # Define test and trial spaces
  VK     = TestFESpace(Ω  , reffeᵤ; conformity=:L2)
  V∂K    = TestFESpace(Γ  , reffeᵤ; conformity=:L2)
  UK     = TrialFESpace(VK)
  U∂K    = TrialFESpace(V∂K)
  VK_V∂K = MultiFieldFESpace([VK,V∂K])
  UK_U∂K = MultiFieldFESpace([UK,U∂K])

  degree = 2*order+1
  dΩ     = Measure(Ω,degree)
  dΓ     = Measure(Γ,degree)
  nK     = get_cell_normal_vector(∂K)
  d∂K    = Measure(∂K,degree)

  reconstruction_op=setup_reconstruction_operator(model, order, dΩ, d∂K)
  uhK_uh∂K =get_trial_fe_basis(UK_U∂K)
  reconstruction_op_image_span=reconstruction_op(uhK_uh∂K)
  uh_dofs=zeros(num_free_dofs(UK_U∂K))
  uh_dofs[2]=1.0
  uhK_uh∂K=FEFunction(UK_U∂K,uh_dofs)
  projected_uhK_uh∂K=reconstruction_op(uhK_uh∂K)
  # (uK,u∂K) = get_trial_fe_basis(UK_U∂K)
  # v_c,v_nc = get_fe_basis(VKR_DS_DECOMP)
  # u_nzm,u_zm = get_trial_fe_basis(UKR_DS_DECOMP)
  # dc=∫((∇(v_nc)⋅nK)*u∂K)d∂K
  # dc=∫(∇(v_nc)⋅∇(u_zm))dΩ + ∫(∇(v_nc)⋅∇(u_nzm))dΩ + ∫(v_c*u_nzm)dΩ
  # dc=∫(v_c*u_zm)dΩ

  reduction_op=setup_reduction_operator(UK_U∂K,VK_V∂K,dΩ,d∂K)
  x=reduction_op(get_trial_fe_basis(UK))
  uhk=FEFunction(UK,rand(num_free_dofs(UK)))
  y=reduction_op(uhk)

  # Check that the result of applying the reconstruction operator to the
  # the result of applyting the reduction operator to P_K^{k+1} is P_K^{k+1}
  # itself
  reffeᵤ     = ReferenceFE(lagrangian,Float64,order+1)
  VH1        = TestFESpace(Ω, reffeᵤ; conformity=:H1)
  UH1        = TrialFESpace(VH1)
  free_dofs  = rand(num_free_dofs(UH1))
  uh         = FEFunction(UH1, free_dofs)
  uh_reduced = reduction_op(uh)

  uh_reconstructed = reconstruction_op(uh_reduced)
  eh = uh_reconstructed-uh
  @test sum(∫(eh*eh)dΩ) < 1.0e-12
end
