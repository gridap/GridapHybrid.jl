module CellBoundaryTests
   using Test
   using Gridap
   using ExploringGridapHybridization
   domain = (0,4,0,4)
   partition = (2,2)
   model = CartesianDiscreteModel(domain,partition)
   D=2
   model_Γ = Gridap.Geometry.BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))
   ∂Tbis=CellBoundaryBis(model)
   m=Gridap.Geometry.get_cell_ref_map(∂Tbis)
   xbis,wbis=quadrature_points_and_weights(∂Tbis,2)
   mx=lazy_map(evaluate,m,xbis)

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

   @time vh_∂Tbis   = restrict_to_cell_boundary(∂Tbis,vh)
   @time vh_∂Tbis_x = lazy_map(evaluate,vh_∂Tbis,xbis)

   ∂T=CellBoundary(model)
   x,w=quadrature_points_and_weights(∂T,2)
   @time vh_∂T   = restrict_to_cell_boundary(∂T,vh)
   @time vh_∂T_x = lazy_map(evaluate,vh_∂T,x)
   @test all(vh_∂Tbis_x == vh_∂T_x)

end
