module CellBoundaryTests
   using Test
   using Gridap
   using ExploringGridapHybridization
   domain = (0,4,0,4)
   partition = (2,2)
   model = CartesianDiscreteModel(domain,partition)
   D=2
   model_Γ = Gridap.Geometry.BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))
   ∂T=CellBoundary(model)
   m=Gridap.Geometry.get_cell_ref_map(∂T)
   x,w=quadrature_points_and_weights(∂T,2)
   mx=lazy_map(evaluate,m,x)

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

   @time vh_∂T   = restrict_to_cell_boundary(∂T,vh)
   @time vh_∂T_x = lazy_map(evaluate,vh_∂T,x)

   @time uh_∂T   = restrict_to_cell_boundary(∂T,uh)
   @time uh_∂T_x = lazy_map(evaluate,uh_∂T,x)

   @time mh_∂T   = restrict_to_cell_boundary(∂T,mh)
   @time mh_∂T_x = lazy_map(evaluate,mh_∂T,x)

   @time lh_∂T   = restrict_to_cell_boundary(∂T,lh)
   @time lh_∂T_x = lazy_map(evaluate,lh_∂T,x)

   @time n=get_cell_normal_vector(∂T)
   @time nx=lazy_map(evaluate,n,x)

   @time nown=get_cell_owner_normal_vector(∂T)
   @time nownx=lazy_map(evaluate,nown,x)

   @time map=get_cell_map(∂T)
   @time mapx=lazy_map(evaluate,map,x)

   @time mapgrad=lazy_map(∇,map)
   @time mapgradx=lazy_map(evaluate,mapgrad,x)


   #print_op_tree(map)
   #mapf  = get_cell_map(∂T.btrian)
   #fmapf = Gridap.Geometry.get_cell_ref_map(∂T.btrian)
   #print_op_tree(fmapf)
   #print_op_tree(lazy_map(Broadcasting(∇),fmapf))
   #print_op_tree(lazy_map(∇,fmapf))
   #print_op_tree(lazy_map(∇,fmapf))
   #x=lazy_map(Broadcasting(∇),map)
end
