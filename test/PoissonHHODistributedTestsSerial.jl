module PoissonHHODistributedTests
  using Gridap
  using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.Arrays, Gridap.CellData, Gridap.FESpaces
  using GridapHybrid
  using GridapDistributed
  using PartitionedArrays
  using Plots

  function setup_reconstruction_operator(model, order, dΩ, d∂K)
    nK        = get_cell_normal_vector(d∂K.quad.trian)
    refferecᵤ = ReferenceFE(orthogonal_basis, Float64, order+1)
    reffe_c   = ReferenceFE(monomial_basis  , Float64, order+1; subspace=:OnlyConstant)

    Ω = dΩ.quad.trian

    VKR     = TestFESpace(Ω  , refferecᵤ; conformity=:L2)
    UKR     = TrialFESpace(VKR)
    VKR_C   = TestFESpace(Ω, reffe_c ; conformity=:L2, vector_type=Vector{Float64})
    UKR_C   = TrialFESpace(VKR_C)

    V = MultiFieldFESpace([VKR,VKR_C])
    U = MultiFieldFESpace([UKR,UKR_C])

    m( (u,u_c), (v,v_c) ) = ∫(∇(v)⋅∇(u))dΩ + ∫(v_c*u)dΩ + ∫(v*u_c)dΩ
    n( (uK,u∂K), (v,v_c) ) = ∫(∇(v)⋅∇(uK))dΩ + ∫(v_c*uK)dΩ + ∫((∇(v)⋅nK)*u∂K)d∂K - ∫((∇(v)⋅nK)*uK)d∂K 

    ReconstructionFEOperator((m,n), U, V)
  end

  function setup_projection_operator(UK_U∂K,VK_V∂K,R,dΩ,d∂K)
    m( (uK,u∂K), (vK,v∂K) ) = ∫(vK*uK)dΩ + ∫(v∂K*u∂K)d∂K
    function n(uK_u∂K, (vK,v∂K))
      urK_ur∂K = R(uK_u∂K)
      urK,ur∂K = urK_ur∂K
      uK,u∂K   = uK_u∂K
      # ∫(vK*(urK-uK))dΩ+∫(vK*ur∂K)dΩ -                 # bulk projection terms
      #    ∫(v∂K*urK)d∂K+∫(v∂K*u∂K)d∂K-∫(v∂K*ur∂K)d∂K   # skeleton projection terms
      ∫(vK*(urK-uK))dΩ+∫(vK*ur∂K)dΩ +                   # bulk projection terms
        ∫(v∂K*urK)d∂K-∫(v∂K*u∂K)d∂K+∫(v∂K*ur∂K)d∂K     # skeleton projection terms
    end
    function n(uK_u∂K::Gridap.MultiField.MultiFieldFEFunction, (vK,v∂K))
      uK,u∂K = uK_u∂K
      urK = R(uK_u∂K)
      ∫(vK*(urK-uK))dΩ +      # bulk projection terms
        ∫(v∂K*(urK-u∂K))d∂K   # skeleton projection terms
    end
    GridapHybrid.ProjectionFEOperator((m,n),UK_U∂K,VK_V∂K,R)
  end

  # Compute length of diagonals in the reference domain
  function foo(m)
  p0 = evaluate(m,Point(0.0,0.0))
  p1 = evaluate(m,Point(1.0,1.0))
  norm(p1-p0)
  # To Do: incorporate the other diagonal    
  #   p2 = evaluate(m,Point(1.0,0.0))
  #   p3 = evaluate(m,Point(0.0,1.0))
  #   maximum( norm(p1-p0), norm(p2-p3) )
  end

  function setup_reconstruction_operator(model::GridapDistributed.GenericDistributedDiscreteModel, 
    order, dΩ::GridapDistributed.DistributedMeasure, d∂K::GridapDistributed.DistributedMeasure)
    R = map(local_views(model), local_views(dΩ), local_views(d∂K)) do m, dΩi, d∂Ki
        setup_reconstruction_operator(m, order, dΩi, d∂Ki)
    end
    return R
  end  

  function setup_projection_operator(trial::GridapDistributed.DistributedMultiFieldFESpace, 
    test::GridapDistributed.DistributedMultiFieldFESpace, R::AbstractArray{<:ReconstructionFEOperator},
    dΩ::GridapDistributed.DistributedMeasure, d∂K::GridapDistributed.DistributedMeasure)
    P = map(local_views(trial), local_views(test), local_views(R), local_views(dΩ), local_views(d∂K)) do triali, testi, Ri, dΩi, d∂Ki
        setup_projection_operator(triali, testi, Ri, dΩi, d∂Ki)
    end
    return P
  end  

  function r(u,v,R,dΩ)
    uK_u∂K=R(u)
    vK_v∂K=R(v)
    uK,u∂K = uK_u∂K
    vK,v∂K = vK_v∂K 
    ∫(∇(vK)⋅∇(uK))dΩ + ∫(∇(vK)⋅∇(u∂K))dΩ +
       ∫(∇(v∂K)⋅∇(uK))dΩ + ∫(∇(v∂K)⋅∇(u∂K))dΩ
  end

  function s(u,v,P,d∂K)
    # Currently, we cannot use this CellField in the integrand of the skeleton integrals below.
    # we think we need to develop a specific version for _transform_face_to_cell_lface_expanded_array
    # we will be using h_T_1 in the meantime
    
    # Ω = get_triangulation(u)
    # cmaps = collect1d(get_cell_map(Ω))
    # h_array = lazy_map(foo,cmaps)
    # # h_T = CellField(h_array,Ω,PhysicalDomain())
    # h_T = CellField(h_array,Ω)
    # h_T_1 = 1.0/h_T

    lΩ = get_triangulation(u)
    dlΩ = Measure(lΩ, 3)
    h_T=CellField(get_array(∫(1)dlΩ),lΩ)
    h_T_1=1.0/h_T    
 
    uK_u∂K_ΠK,uK_u∂K_Π∂K=P(u)
    vK_v∂K_ΠK,vK_v∂K_Π∂K=P(v)

    uK_ΠK  , u∂K_ΠK  = uK_u∂K_ΠK
    uK_Π∂K , u∂K_Π∂K = uK_u∂K_Π∂K
    
    vK_ΠK  , v∂K_ΠK  = vK_v∂K_ΠK
    vK_Π∂K , v∂K_Π∂K = vK_v∂K_Π∂K

    vK_Π∂K_vK_ΠK=vK_Π∂K-vK_ΠK
    v∂K_Π∂K_v∂K_ΠK=v∂K_Π∂K-v∂K_ΠK
    uK_Π∂K_uK_ΠK=uK_Π∂K-uK_ΠK
    u∂K_Π∂K_u∂K_ΠK=u∂K_Π∂K-u∂K_ΠK

    ∫(h_T_1*(vK_Π∂K_vK_ΠK)*(uK_Π∂K_uK_ΠK))d∂K + 
        ∫(h_T_1*(v∂K_Π∂K_v∂K_ΠK)*(u∂K_Π∂K_u∂K_ΠK))d∂K +
         ∫(h_T_1*(vK_Π∂K_vK_ΠK)*(u∂K_Π∂K_u∂K_ΠK))d∂K + 
          ∫(h_T_1*(v∂K_Π∂K_v∂K_ΠK)*(uK_Π∂K_uK_ΠK))d∂K
  end

  function r(u::GridapDistributed.DistributedMultiFieldCellField, 
             v::GridapDistributed.DistributedMultiFieldCellField,
             R::AbstractArray{<:ReconstructionFEOperator},
             dΩ::GridapDistributed.DistributedMeasure)
    a_r = map(local_views(u), local_views(v), local_views(R), local_views(dΩ)) do ui, vi, Ri, dΩi
                  a_r = r(ui,vi,Ri,dΩi)
        return a_r
    end
  end

  function s(u::GridapDistributed.DistributedMultiFieldCellField, 
             v::GridapDistributed.DistributedMultiFieldCellField,
             P::AbstractArray{<:ProjectionFEOperator}, 
             d∂K::GridapDistributed.DistributedMeasure)
    a_s = map(local_views(u), local_views(v), local_views(P), local_views(d∂K)) do ui, vi, Pi, d∂Ki
        a_s = s(ui,vi,Pi,d∂Ki)
        return a_s
    end
  end

  function solve_hho(domain, mesh_partition, np, ranks, order, u)  
    model = CartesianDiscreteModel(ranks, np, domain, mesh_partition)

    D = num_cell_dims(model)
    Ω = Triangulation(ReferenceFE{D},model)
    Γ = Triangulation(ReferenceFE{D-1},model)
    ∂K= GridapHybrid.Skeleton(model)

    # Local Triangulations to be used for local problems
    ∂K_loc = GridapDistributed.DistributedTriangulation(
      map(t -> t.parent,local_views(∂K)), ∂K.model)
    Ω_loc = GridapDistributed.add_ghost_cells(Ω)

    reffe = ReferenceFE(lagrangian, Float64, order; space=:P)
    VK    = TestFESpace(Ω  , reffe; conformity=:L2)
    V∂K   = TestFESpace(Γ  , reffe; conformity=:L2, dirichlet_tags=collect(5:8))
    UK = TrialFESpace(VK)
    U∂K= TrialFESpace(V∂K,u)
    V = MultiFieldFESpace([VK,V∂K])
    U = MultiFieldFESpace([UK,U∂K])

    degree = 2*(order+1)
    dΩ = Measure(Ω, degree)
    d∂K= Measure(∂K, degree)
    d∂K_loc = Measure(∂K_loc, degree)
    dΩ_loc = Measure(Ω_loc, degree)
    # Ω_loc = get_triangulation(VK)  # alternative way to define Ω_loc which ensures that the object id is same as that of VK.
    # ndofs= num_free_dofs(U)

    R = setup_reconstruction_operator(model, order, dΩ_loc, d∂K_loc)
    P = setup_projection_operator(U, V, R, dΩ_loc, d∂K_loc)
    
    a(u,v,dΩ,d∂K) = map(+, local_views(r(u,v,R,dΩ)), local_views(s(u,v,P,d∂K)))
    l((vK,),dΩ) = ∫(vK*f)dΩ

    weakform = (u,v)->(a(u,v,dΩ,d∂K),l(v,dΩ))
    op=HybridAffineFEOperator(weakform, U, V, [1], [2])

    xh = solve(op);
    uhK,uh∂K = xh
    e = uhK-u

    return sqrt(sum(∫(e⋅e)dΩ))
  end

  function convg_test(domain,ns,np,order,u)
    
    ranks = with_debug() do distribute
      distribute(LinearIndices((prod(np),)))
    end
    
    el2 = Float64[]
    hs = Float64[]
    for n in ns
      l2 = solve_hho(domain,(n,n),np,ranks,order,u)
      println(l2)
      push!(el2,l2)
      h = 1/n
      push!(hs,h)
    end
    println(el2)
    println(hs)
    el2, hs
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x), 1), x) \ y
    linreg[2]
  end


  u(x) = sin(2π*x[1])*sin(2π*x[2])
  f(x) = -Δ(u)(x)

  domain = (0,1,0,1)

  np = (2,2)
  #ns = [np[1]*2^i for i in 2:5]
  ns = [2^i for i in 2:3]
  # mesh_partition = (2,2)

  order = 1
  el, hs = convg_test(domain,ns,np,order,u)
  println("Slope L2-norm u: $(slope(hs,el))")
  slopekp1=[Float64(ni)^(-(order+1)) for ni in ns]
  display(plot(hs,[el slopekp1],
    xaxis=:log, yaxis=:log,
    label=["L2u (measured)" "slope k+1"],
    shape=:auto,
    xlabel="h",ylabel="L2 error",legend=:bottomright))
  
end 
