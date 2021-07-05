module CellBoundaryTests
   using Test
   using Gridap
   using ExploringGridapHybridization
   domain = (0,4,0,4)
   partition = (2,2)
   model = CartesianDiscreteModel(domain,partition)
   D=2
   model_Γ = Gridap.Geometry.BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))
   ∂Topt=CellBoundaryOpt(model)
   m=Gridap.Geometry.get_cell_ref_map(∂Topt)
   xopt,wopt=quadrature_points_and_weights(∂Topt,2)
   mx=lazy_map(evaluate,m,xopt)

   # Define test FESpaces
   order=0
   reffeᵤ = ReferenceFE(raviart_thomas,Float64,order)
   reffeₚ = ReferenceFE(lagrangian,Float64,order)
   reffeₗ = ReferenceFE(lagrangian,Float64,order)
   V = TestFESpace(model  , reffeᵤ; conformity=:L2)
   Q = TestFESpace(model  , reffeₚ; conformity=:L2)
   M = TestFESpace(model_Γ, reffeₗ; conformity=:L2)
   Y = MultiFieldFESpace([V,Q,M])

   # Create trial spaces
   U = TrialFESpace(V)
   P = TrialFESpace(Q)
   L = TrialFESpace(M)
   X = MultiFieldFESpace([U, P, L])

   yh = get_fe_basis(Y)
   vh,qh,mh = yh

   xh = get_trial_fe_basis(X)
   uh,ph,lh = xh

   ∂T=CellBoundary(model)
   x,w=quadrature_points_and_weights(∂T,2)

   @time vh_∂Topt   = restrict_to_cell_boundary(∂Topt,vh)
   @time vh_∂Topt_x = lazy_map(evaluate,vh_∂Topt,xopt)

   @time vh_∂T   = restrict_to_cell_boundary(∂T,vh)
   @time vh_∂T_x = lazy_map(evaluate,vh_∂T,x)
   @test all(vh_∂Topt_x == vh_∂T_x)

   @time uh_∂Topt   = restrict_to_cell_boundary(∂Topt,uh)
   @time uh_∂Topt_x = lazy_map(evaluate,uh_∂Topt,xopt)

   @time uh_∂T   = restrict_to_cell_boundary(∂T,uh)
   @time uh_∂T_x = lazy_map(evaluate,uh_∂T,x)
   @test all(uh_∂Topt_x == uh_∂T_x)

   @time mh_∂Topt   = restrict_to_cell_boundary(∂Topt,mh)
   @time mh_∂Topt_x = lazy_map(evaluate,mh_∂Topt,x)

   @time mh_∂T   = restrict_to_cell_boundary(∂T,mh)
   @time mh_∂T_x = lazy_map(evaluate,mh_∂T,x)
   @test all(mh_∂Topt_x == mh_∂T_x)

   @time lh_∂Topt   = restrict_to_cell_boundary(∂Topt,lh)
   @time lh_∂Topt_x = lazy_map(evaluate,lh_∂Topt,x)

   @time lh_∂T   = restrict_to_cell_boundary(∂T,lh)
   @time lh_∂T_x = lazy_map(evaluate,lh_∂T,x)
   @test all(lh_∂Topt_x == lh_∂T_x)

   # nf=Gridap.Geometry.get_facet_normal(∂Topt.btrian)
   # dΓ=Measure(∂Topt.btrian,2)
   # nfx=lazy_map(evaluate,nf,Gridap.CellData.get_data(get_cell_points(dΓ.quad)))

   @time n=get_cell_normal_vector(∂T)
   @time nx=lazy_map(evaluate,n,x)

   @time nopt=get_cell_normal_vector(∂Topt)
   @time noptx=lazy_map(evaluate,n,xopt)

   @test all(nx == noptx)
end
