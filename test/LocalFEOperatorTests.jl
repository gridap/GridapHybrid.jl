module LocalFEOperatorTests
  using Gridap
  using GridapHybrid
  using Test
  using LinearAlgebra

  model=CartesianDiscreteModel((0,1,0,1),(2,2))

  # Geometry
  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D},model)
  Γ = Triangulation(ReferenceFE{D-1},model)
  ∂K = GridapHybrid.Skeleton(model)

  # Reference FEs
  order     = 1
  T         = Float64
  reffe_nzm = ReferenceFE(orthogonal_basis, T, order+1; subspace=:NonZeroMean)
  reffe_zm  = ReferenceFE(orthogonal_basis, T, order+1; subspace=:ZeroMean)
  reffe_c   = ReferenceFE(monomial_basis, T, order+1; subspace=:OnlyConstant)
  reffe_nc  = ReferenceFE(monomial_basis, T, order+1; subspace=:ExcludeConstant)

  UKR_NZM = TrialFESpace(TestFESpace(Ω, reffe_nzm; conformity=:L2))
  UKR_ZM  = TrialFESpace(TestFESpace(Ω, reffe_zm; conformity=:L2))
  VKR_C   = TestFESpace(Ω, reffe_c ; conformity=:L2, vector_type=Vector{Float64})
  VKR_NC  = TestFESpace(Ω, reffe_nc; conformity=:L2, vector_type=Vector{Float64})

  VKR_DS_DECOMP = MultiFieldFESpace([VKR_C,VKR_NC])
  UKR_DS_DECOMP = MultiFieldFESpace([UKR_NZM,UKR_ZM])

  refferecᵤ = ReferenceFE(orthogonal_basis,Float64,order+1)
  reffeᵤ    = ReferenceFE(lagrangian,Float64,order  ;space=:P)

  VKR  = TestFESpace(Ω  , refferecᵤ; conformity=:L2)
  VK   = TestFESpace(Ω  , reffeᵤ   ; conformity=:L2)
  V∂K  = TestFESpace(Γ  , reffeᵤ   ; conformity=:L2)

  VK_V∂K = MultiFieldFESpace([VK,V∂K])

  # Define trial FEspaces
  UKR = TrialFESpace(VKR)
  UK  = TrialFESpace(VK)
  U∂K = TrialFESpace(V∂K)

  UK_U∂K   = MultiFieldFESpace([UK,U∂K])

  degree = 2*order+1
  dΩ     = Measure(Ω,degree)
  nK     = get_cell_normal_vector(∂K)
  d∂K    = Measure(∂K,degree)


  # (uK,u∂K) = get_trial_fe_basis(UK_U∂K)
  # v_c,v_nc = get_fe_basis(VKR_DS_DECOMP)
  # u_nzm,u_zm = get_trial_fe_basis(UKR_DS_DECOMP)
  # dc=∫((∇(v_nc)⋅nK)*u∂K)d∂K
  # dc=∫(∇(v_nc)⋅∇(u_zm))dΩ + ∫(∇(v_nc)⋅∇(u_nzm))dΩ + ∫(v_c*u_nzm)dΩ
  # dc=∫(v_c*u_zm)dΩ

  m( (u_nzm,u_zm), (v_c,v_nc) )= ∫(∇(v_nc)⋅∇(u_zm))dΩ +
                                 ∫(∇(v_nc)⋅∇(u_nzm))dΩ +
                                 ∫(v_c*u_nzm)dΩ
  n( (uK,u∂K), (v_c,v_nc)     )= ∫(-Δ(v_nc)*uK)dΩ + ∫(v_c*uK)dΩ + ∫((∇(v_nc)⋅nK)*u∂K)d∂K


  op=LocalFEOperator((m,n),UKR,VKR;
                     trial_space_ds_decomp=UKR_DS_DECOMP,
                     test_space_ds_decomp=VKR_DS_DECOMP)
  UKR_basis=op(get_trial_fe_basis(UK_U∂K))
  uh_dofs=zeros(num_free_dofs(UK_U∂K))
  uh_dofs[2]=1.0
  puh=UKR_basis=op(FEFunction(UK_U∂K,uh_dofs))
end
