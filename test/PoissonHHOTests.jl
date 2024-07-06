# module PoissonHHOTests
  using Gridap
  using GridapHybrid
  using Test
  using Plots

  function setup_reconstruction_operator(model, order, dΩ, d∂K, VK_V∂K)
    nK        = get_cell_normal_vector(d∂K.quad.trian)
    refferecᵤ = ReferenceFE(orthogonal_basis, Float64, order+1)
    
    # reffe_nzm = ReferenceFE(orthogonal_basis, Float64, order+1; subspace=:NonZeroMean)
    # reffe_zm  = ReferenceFE(orthogonal_basis, Float64, order+1; subspace=:ZeroMean)
    reffe_c   = ReferenceFE(monomial_basis  , Float64, order+1; subspace=:OnlyConstant)
    # reffe_nc  = ReferenceFE(monomial_basis  , Float64, order+1; subspace=:ExcludeConstant)

    Ω = dΩ.quad.trian

    VKR     = TestFESpace(Ω  , refferecᵤ; conformity=:L2)
    UKR     = TrialFESpace(VKR)
    # UKR_NZM = TrialFESpace(TestFESpace(Ω, reffe_nzm; conformity=:L2))
    # UKR_ZM  = TrialFESpace(TestFESpace(Ω, reffe_zm; conformity=:L2))
    VKR_C   = TestFESpace(Ω, reffe_c ; conformity=:L2, vector_type=Vector{Float64})
    UKR_C   = TrialFESpace(VKR_C)
    # VKR_NC  = TestFESpace(Ω, reffe_nc; conformity=:L2, vector_type=Vector{Float64})

    # VKR_DS_DECOMP = MultiFieldFESpace([VKR_C,VKR_NC])
    # UKR_DS_DECOMP = MultiFieldFESpace([UKR_NZM,UKR_ZM])

    V = MultiFieldFESpace([VKR,VKR_C])
    U = MultiFieldFESpace([UKR,UKR_C])

    # m( (u_nzm,u_zm), (v_c,v_nc) ) = ∫(∇(v_nc)⋅∇(u_zm))dΩ + ∫(∇(v_nc)⋅∇(u_nzm))dΩ +
    #                                ∫(v_c*u_nzm)dΩ
    # n( (uK,u∂K), (v_c,v_nc)     ) = ∫(-Δ(v_nc)*uK)dΩ + ∫(v_c*uK)dΩ + ∫((∇(v_nc)⋅nK)*u∂K)d∂K

    m( (u,u_c), (v,v_c) ) = ∫(∇(v)⋅∇(u))dΩ + ∫(v_c*u)dΩ + ∫(v*u_c)dΩ
    n( (uK,u∂K), (v,v_c) ) = ∫(∇(v)⋅∇(uK))dΩ + ∫(v_c*uK)dΩ + ∫((∇(v)⋅nK)*u∂K)d∂K - ∫((∇(v)⋅nK)*uK)d∂K 

    ReconstructionFEOperator((m,n), U, V)

    # LocalFEOperator((m,n),UKR,VKR;
    #                 trial_space_ds_decomp=UKR_DS_DECOMP,
    #                 test_space_ds_decomp=VKR_DS_DECOMP)
  end

  function setup_projection_operator(UK_U∂K,VK_V∂K,R,dΩ,d∂K)
    m((uK,u∂K)  , (vK,v∂K)) = ∫(vK*uK)dΩ + ∫(v∂K*u∂K)d∂K
    function n(uK_u∂K, (vK,v∂K))
      urK_ur∂K = R(uK_u∂K)
      urK,ur∂K = urK_ur∂K
      uK,u∂K   = uK_u∂K
      ∫(vK*(urK-uK))dΩ-∫(vK*ur∂K)dΩ -                 # bulk projection terms
         ∫(v∂K*urK)d∂K+∫(v∂K*u∂K)d∂K-∫(v∂K*ur∂K)d∂K   # skeleton projection terms
    end
    ProjectionFEOperator((m,n),UK_U∂K,VK_V∂K)
   end

   p = 1
   u(x) = x[1]^p+x[2]^p                         # Ex 1
  #  u(x) = x[1]*(x[1]-1)^p*x[2]*(x[2]-1)^p         # Ex 2
   f(x)=-Δ(u)(x)

   #u(x)=x[1]+x[2]
   #f(x)=-Δ(u)(x)

  # function solve_hho(cells,order)
      cells=(1,1)
      order=0
      partition = (0,1,0,1)
      model = CartesianDiscreteModel(partition, cells)
      D  = num_cell_dims(model)
      Ω  = Triangulation(ReferenceFE{D},model)
      Γ  = Triangulation(ReferenceFE{D-1},model)
      ∂K = GridapHybrid.Skeleton(model)

      reffeᵤ    = ReferenceFE(lagrangian,Float64,order  ;space=:P)

      VK     = TestFESpace(Ω  , reffeᵤ; conformity=:L2)
      V∂K    = TestFESpace(Γ  , reffeᵤ; conformity=:L2,dirichlet_tags=collect(5:8))
      UK     = TrialFESpace(VK)
      U∂K    = TrialFESpace(V∂K,u)
      VK_V∂K = MultiFieldFESpace([VK,V∂K])
      UK_U∂K = MultiFieldFESpace([UK,U∂K])

      degree = 2*order+1
      dΩ     = Measure(Ω,degree)
      dΓ     = Measure(Γ,degree)
      d∂K    = Measure(∂K,degree)

      R=setup_reconstruction_operator(model, order, dΩ, d∂K, VK_V∂K)
      projection_op=setup_projection_operator(UK_U∂K,VK_V∂K,R,dΩ,d∂K)

      function r(u,v)
        uK_u∂K=R(u)
        vK_v∂K=R(v)
        uK,u∂K = uK_u∂K
        vK,v∂K = vK_v∂K 
        ∫(∇(vK)⋅∇(uK))dΩ + ∫(∇(vK)⋅∇(u∂K))dΩ +
           ∫(∇(v∂K)⋅∇(uK))dΩ + ∫(∇(v∂K)⋅∇(u∂K))dΩ
      end

      function s(u,v)
        h_T=CellField(get_array(∫(1)dΩ),Ω)
        h_T_1=1.0/h_T
        h_T_2=1.0/(h_T*h_T)

        # Currently, I cannot use this CellField in the integrand of the skeleton integrals below.
        # I think we need to develop a specific version for _transform_face_to_cell_lface_expanded_array
        # I will be using h_T_1 in the meantime
        h_F=CellField(get_array(∫(1)dΓ),Γ)
        h_F=1.0/h_F

        uK_u∂K_ΠK,uK_u∂K_Π∂K=projection_op(u)
        vK_v∂K_ΠK,vK_v∂K_Π∂K=projection_op(v)

        uK_ΠK  , u∂K_ΠK  = uK_u∂K_ΠK
        uK_Π∂K , u∂K_Π∂K = uK_u∂K_Π∂K
        
        vK_ΠK  , v∂K_ΠK  = vK_v∂K_ΠK
        vK_Π∂K , v∂K_Π∂K = vK_v∂K_Π∂K

        ∫(h_T_1*(vK_Π∂K-vK_ΠK)*(uK_Π∂K-uK_ΠK))d∂K + 
           ∫(h_T_1*(v∂K_Π∂K-v∂K_ΠK)*(u∂K_Π∂K-u∂K_ΠK))d∂K +
            ∫(h_T_1*(vK_Π∂K-vK_ΠK)*(u∂K_Π∂K-u∂K_ΠK))d∂K + 
             ∫(h_T_1*(v∂K_Π∂K-v∂K_ΠK)*(uK_Π∂K-uK_ΠK))d∂K
      end

      a(u,v)=r(u,v)+s(u,v)
      l((vK,))=∫(vK*f)dΩ


      op=HybridAffineFEOperator((u,v)->(a(u,v),l(v)), UK_U∂K, VK_V∂K, [1], [2])
      xh=solve(op)

      uhK,uh∂K=xh

      e = u -uhK
      # @test sqrt(sum(∫(e⋅e)dΩ)) < 1.0e-12
      return sqrt(sum(∫(e⋅e)dΩ))
  end

  function conv_test(ns,order)
    el2 = Float64[]
    hs = Float64[]
    for n in ns
      l2 = solve_hho((n,n),order)
      println(l2)
      h = 1.0/n
      push!(el2,l2)
      push!(hs,h)
    end
    println(el2)
    el2, hs
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x), 1), x) \ y
    linreg[2]
  end

  solve_hho((2,2),0)

  # ns=[8,16,32,64,128]
  ns=[8,16,32,64]
  order=0
  el, hs = conv_test(ns,order)
  println("Slope L2-norm u: $(slope(hs,el))")
  slopek  =[Float64(ni)^(-(order)) for ni in ns]
  slopekp1=[Float64(ni)^(-(order+1)) for ni in ns]
  slopekp2=[Float64(ni)^(-(order+2)) for ni in ns]
  display(plot(hs,[el slopek slopekp1 slopekp2],
    xaxis=:log, yaxis=:log,
    label=["L2u (measured)" "slope k" "slope k+1" "slope k+2"],
    shape=:auto,
    xlabel="h",ylabel="L2 error",legend=:bottomright))
end
