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

  function setup_difference_operator(UK_U∂K,VK_V∂K,R,dΩ,d∂K)
    m((uK,u∂K)  , (vK,v∂K)) = ∫(vK*uK)dΩ + ∫(v∂K*u∂K)d∂K
    function n(uK_u∂K, (vK,v∂K))
      uK_u∂K_rec     = R(uK_u∂K)
      uK,u∂K         = uK_u∂K
      uK_rec,u∂K_rec = uK_u∂K_rec
      ∫(vK *uK_rec)dΩ  + ∫(vK*u∂K_rec)dΩ   - ∫(vK *uK)dΩ +
      ∫(v∂K*uK_rec)d∂K + ∫(v∂K*u∂K_rec)d∂K - ∫(v∂K*u∂K)d∂K
    end
    LocalFEOperator((m,n),UK_U∂K,VK_V∂K; field_type_at_common_faces=MultiValued())
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

  R=setup_reconstruction_operator(model, order, dΩ, d∂K)
  uhK_uh∂K =get_trial_fe_basis(UK_U∂K)
  reconstruction_op_image_span=R(uhK_uh∂K)
  uh_dofs=zeros(num_free_dofs(UK_U∂K))
  uh_dofs[2]=1.0
  uhK_uh∂K=FEFunction(UK_U∂K,uh_dofs)
  projected_uhK_uh∂K=R(uhK_uh∂K)
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
  # the result of applying the reduction operator to P_K^{k+1} is P_K^{k+1}
  # itself
  reffeᵤ     = ReferenceFE(lagrangian,Float64,order+1)
  VH1        = TestFESpace(Ω, reffeᵤ; conformity=:H1)
  UH1        = TrialFESpace(VH1)
  free_dofs  = rand(num_free_dofs(UH1))
  uh         = FEFunction(UH1, free_dofs)
  uh_reduced = reduction_op(uh)

  uh_reconstructed = R(uh_reduced)
  eh = uh_reconstructed-uh
  @test sum(∫(eh*eh)dΩ) < 1.0e-12

  # Check that the mean value of the reconstructed cell FE space functions
  # match the mean value of the original cell FE space functions
  uhK_uh∂K =get_trial_fe_basis(UK_U∂K)
  vhK_vh∂K =get_fe_basis(VK_V∂K)
  vK,v∂K=vhK_vh∂K
  uK,u∂K=uhK_uh∂K

  R_vhK_vh∂K=R(vhK_vh∂K)
  R_vhK,_=R_vhK_vh∂K

  dc1=∫(vK)dΩ
  dc2=∫(R_vhK)dΩ

  @test all(get_array(dc1) .≈ get_array(dc2))

  diff_op=setup_difference_operator(UK_U∂K,VK_V∂K,R,dΩ,d∂K)
  v=get_fe_basis(VK_V∂K)
  ub=get_trial_fe_basis(UK_U∂K)

  ub_rec=R(ub)
  ub_rec1,ub_rec2=ub_rec
  xΩ=Gridap.CellData.get_cell_points(dΩ.quad)
  x∂K=Gridap.CellData.get_cell_points(d∂K.quad)

  v_rec=R(v)
  v_rec1,v_rec2=v_rec
  v_rec1_d∂K=Gridap.CellData.change_domain(v_rec1,d∂K.quad.trian,Gridap.CellData.ReferenceDomain())
  v_rec1_d∂K(x∂K)[1][4][1][1]
  v_rec2_d∂K=Gridap.CellData.change_domain(v_rec2,d∂K.quad.trian,Gridap.CellData.ReferenceDomain())
  v_rec2_d∂K(x∂K)[1][3][2][1]

  uK_u∂K_rec     = R(ub)
  uK,u∂K         = ub
  uK_rec,u∂K_rec = uK_u∂K_rec
  dc=∫(vK *uK_rec)dΩ  + ∫(vK*u∂K_rec)dΩ   - ∫(vK *uK)dΩ +
     ∫(v∂K*uK_rec)d∂K + ∫(v∂K*u∂K_rec)d∂K - ∫(v∂K*u∂K)d∂K


  δvK,δv∂K=diff_op(v)
  δuK,δu∂K=diff_op(ub)

  δv∂K_K,δv∂K_∂K=δv∂K
  δu∂K_K,δu∂K_∂K=δu∂K

  δv∂K_K(x∂K)[1][1][1]
  δv∂K_K(x∂K)[1][2][1]
  δv∂K_K(x∂K)[1][3][1]
  δv∂K_K(x∂K)[1][4][1]
  δv∂K_∂K(x∂K)[1][1][2]
  δv∂K_∂K(x∂K)[1][2][2]
  δv∂K_∂K(x∂K)[1][3][2]
  δv∂K_∂K(x∂K)[1][4][2]

  (δv∂K_K*δu∂K_K)(x∂K)[1][4][1,1]
  (δv∂K_K*δu∂K_∂K)(x∂K)[1][4][1,2]
  (δv∂K_∂K*δu∂K_K)(x∂K)[1][4][2,1]
  (δv∂K_∂K*δu∂K_∂K)(x∂K)[1][4][2,2]

  dc=∫(δv∂K_K*δu∂K_K +
      δv∂K_K*δu∂K_∂K+
       δv∂K_∂K*δu∂K_K+
         δv∂K_∂K*δu∂K_∂K)d∂K

  get_array(dc)[1][1,1]
  get_array(dc)[1][2,1]
  get_array(dc)[1][1,2]
  get_array(dc)[1][2,2]
  δvK_K,δvK_∂K=δvK
  δuK_K,δuK_∂K=δuK

  dc=∫(δvK_K*δuK_K)dΩ+
       ∫(δvK_K*δuK_∂K)dΩ+
         ∫(δvK_∂K*δuK_K)dΩ+
           ∫(δvK_∂K*δuK_∂K)dΩ

  get_array(dc)[1][1,1]
  get_array(dc)[1][1,2]
  get_array(dc)[1][2,1]
  get_array(dc)[1][2,2]


end
