module CellBoundaryTests
   using Test
   using Gridap
   using ExploringGridapHybridization
   domain = (0,4,0,4)
   partition = (2,2)
   model = CartesianDiscreteModel(domain,partition)
   D=2
   model_Γ = Gridap.Geometry.BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))
   ∂Tbis=CellBoundaryOpt(model)
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

   ∂T=CellBoundary(model)
   x,w=quadrature_points_and_weights(∂T,2)

   @time vh_∂Tbis   = restrict_to_cell_boundary(∂Tbis,vh)
   @time vh_∂Tbis_x = lazy_map(evaluate,vh_∂Tbis,xbis)

   @time vh_∂T   = restrict_to_cell_boundary(∂T,vh)
   @time vh_∂T_x = lazy_map(evaluate,vh_∂T,x)
   @test all(vh_∂Tbis_x == vh_∂T_x)

   @time uh_∂Tbis   = restrict_to_cell_boundary(∂Tbis,uh)
   @time uh_∂Tbis_x = lazy_map(evaluate,uh_∂Tbis,xbis)

   @time uh_∂T   = restrict_to_cell_boundary(∂T,uh)
   @time uh_∂T_x = lazy_map(evaluate,uh_∂T,x)
   @test all(uh_∂Tbis_x == uh_∂T_x)

   @time mh_∂Tbis   = restrict_to_cell_boundary(∂Tbis,mh)
   @time mh_∂Tbis_x = lazy_map(evaluate,mh_∂Tbis,x)

   @time mh_∂T   = restrict_to_cell_boundary(∂T,mh)
   @time mh_∂T_x = lazy_map(evaluate,mh_∂T,x)
   @test all(mh_∂Tbis_x == mh_∂T_x)

   @time lh_∂Tbis   = restrict_to_cell_boundary(∂Tbis,lh)
   @time lh_∂Tbis_x = lazy_map(evaluate,lh_∂Tbis,x)

   @time lh_∂T   = restrict_to_cell_boundary(∂T,lh)
   @time lh_∂T_x = lazy_map(evaluate,lh_∂T,x)
   @test all(lh_∂Tbis_x == lh_∂T_x)

   # nf=Gridap.Geometry.get_facet_normal(∂Tbis.btrian)
   # dΓ=Measure(∂Tbis.btrian,2)
   # nfx=lazy_map(evaluate,nf,Gridap.CellData.get_data(get_cell_points(dΓ.quad)))

   @time n=get_cell_normal_vector(∂T)
   @time nx=lazy_map(evaluate,n,x)

   @time nbis=get_cell_normal_vector(∂Tbis)
   @time nbisx=lazy_map(evaluate,n,xbis)

   @test all(nx == nbisx)
end
