module DarcyRTHTests
  using Gridap
  using GridapHybrid
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

  function solve_darcy_rth(model,order)
      D = num_cell_dims(model)
      Ω = Triangulation(ReferenceFE{D},model)
      Γ = Triangulation(ReferenceFE{D-1},model)
      ∂K = GridapHybrid.Skeleton(model)

      reffeᵤ = ReferenceFE(raviart_thomas,Float64,order)
      reffeₚ = ReferenceFE(lagrangian,Float64,order)
      reffeₗ = ReferenceFE(lagrangian,Float64,order)

      # Define test FESpaces
      V = TestFESpace(Ω  , reffeᵤ; conformity=:L2)
      Q = TestFESpace(Ω  , reffeₚ; conformity=:L2)
      M = TestFESpace(Γ,
                      reffeₗ;
                      conformity=:L2,
                      dirichlet_tags=collect(5:8))
      Y = MultiFieldFESpace([V,Q,M])

      # Create trial spaces
      U = TrialFESpace(V)
      P = TrialFESpace(Q)
      L = TrialFESpace(M,p)
      X = MultiFieldFESpace([U, P, L])

      yh = get_fe_basis(Y)
      vh,qh,mh = yh

      xh = get_trial_fe_basis(X)
      uh,ph,lh = xh

      degree = 2*(order+1)
      dΩ     = Measure(Ω,degree)

      n   = get_cell_normal_vector(∂K)
      d∂K = Measure(∂K,degree)

      # dc1=∫( vh⋅uh - (∇⋅vh)*ph + qh*(∇⋅uh) )dΩ
      # dc2=∫((vh⋅n)*lh)d∂K
      # dc3=∫(mh*(uh⋅n))d∂K

      a((uh,ph,lh),(vh,qh,mh)) = ∫( vh⋅uh - (∇⋅vh)*ph + qh*(∇⋅uh) )dΩ +
                                ∫((vh⋅n)*lh)d∂K +
                                ∫(mh*(uh⋅n))d∂K
      l((vh,qh,mh)) = ∫( vh⋅f + qh*(∇⋅u))*dΩ

      op=HybridAffineFEOperator((u,v)->(a(u,v),l(v)), X, Y, [1,2], [3])

      xh=solve(op)

      uh,_=xh
      e = u -uh
      sqrt(sum(∫(e⋅e)dΩ))
  end

  partition = (0,1,0,1)
  cells = (2,2)
  model = CartesianDiscreteModel(partition,cells)
  order = 0
  @test solve_darcy_rth(model,order) < 1.0e-12

  model = simplexify(CartesianDiscreteModel(partition,cells))
  @test_broken solve_darcy_rth(model,order) < 1.0e-12

end
