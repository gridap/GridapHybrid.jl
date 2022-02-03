module MultiFieldLagrangeMultipliersTests
  using Gridap
  using ExploringGridapHybridization
  using Test

  #2D problem
  u(x) = VectorValue(1+x[1],1+x[2])
  Gridap.divergence(::typeof(u)) = (x) -> 2
  p(x) = -3.14
  ∇p(x) = zero(x)
  Gridap.∇(::typeof(p)) = ∇p
  f(x) = u(x) + ∇p(x)
  # Normal component of u(x) on Neumann boundary
  function g(x)
    tol=1.0e-14
    if (abs(x[2])<tol)
      return -x[2]
    elseif (abs(x[2]-1.0)<tol)
      return x[2]
    end
    Gridap.Helpers.@check false
  end

  function Gridap.CellData.get_triangulation(f::Gridap.MultiField.MultiFieldCellField)
    s1 = first(f.single_fields)
    trian = get_triangulation(s1)
    trian
  end

  partition = (0,1,0,1)
  cells = (2,2)
  model = CartesianDiscreteModel(partition,cells)

  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D},model)
  Γ = Triangulation(ReferenceFE{D-1},model)
  ∂K = TempSkeleton(model)

  order=0
  reffeᵤ = ReferenceFE(raviart_thomas,Float64,order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order)
  reffeₗ = ReferenceFE(lagrangian,Float64,order)

  # Define test FESpaces
  V1 = TestFESpace(Ω, reffeᵤ; conformity=:L2)
  V2 = TestFESpace(Ω, reffeᵤ; conformity=:L2)
  Q1 = TestFESpace(Ω, reffeₚ; conformity=:L2)
  Q2 = TestFESpace(Ω, reffeₚ; conformity=:L2)
  M1 = TestFESpace(Γ,
                  reffeₗ;
                  conformity=:L2,
                  dirichlet_tags=collect(5:8))
  M2 = TestFESpace(Γ,
                  reffeₗ;
                  conformity=:L2,
                  dirichlet_tags=collect(5:8))
  Y = MultiFieldFESpace([V1,V2,Q1,Q2,M1,M2])

  # Create trial spaces
  U1 = TrialFESpace(V1)
  U2 = TrialFESpace(V2)
  P1 = TrialFESpace(Q1)
  P2 = TrialFESpace(Q2)
  L1 = TrialFESpace(M1,p)
  L2 = TrialFESpace(M2,p)
  X = MultiFieldFESpace([U1,U2,P1,P2,L1,L2])

  order  = 0
  degree = 2*(order+1)
  dΩ     = Measure(Ω,degree)

  n   = get_cell_normal_vector(∂K)
  d∂K = Measure(∂K,degree)

  ablock((uh,ph,lh),(vh,qh,mh)) = ∫( vh⋅uh - (∇⋅vh)*ph + qh*(∇⋅uh) )dΩ +
                            ∫((vh⋅n)*lh)d∂K +
                            ∫(mh*(uh⋅n))d∂K
  lblock((vh,qh,mh)) = ∫( vh⋅f + qh*(∇⋅u))*dΩ

  a((uh1,uh2,ph1,ph2,lh1,lh2),(vh1,vh2,qh1,qh2,mh1,mh2))=
       ablock((uh1,ph1,lh1),(vh1,qh1,mh1))+ablock((uh2,ph2,lh2),(vh2,qh2,mh2))
  l((vh1,vh2,qh1,qh2,mh1,mh2))=lblock((vh1,qh1,mh1))+lblock((vh2,qh2,mh2))

  op=HybridAffineFEOperator((u,v)->(a(u,v),l(v)), X, Y, collect(1:4), collect(5:6))

  xh=solve(op)

  uh1,uh2,_=xh
  e = (u-uh1)+(u-uh2)
  @test sqrt(sum(∫(e⋅e)dΩ)) < 1.0e-12

end
