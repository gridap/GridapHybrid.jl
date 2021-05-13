# Params
κ = 2.0
κinv = 1/k
c = 3.0
# Manufactured solution
u(x) = x[1] + x[2]
q(x) = -κ*∇(u)(x)
f(x) = ∇⋅(y->c*u(y)+q(y))(x)
# Model
domain = (0,1,0,1)
cells = (3,3)
D = length(cells)
model = CartesianDiscreteModel(domain,cells)
model_facets = DiscreteModel(Polytope{D-1},model)
# FESpaces
k = 5
reffe_w = ReferenceFE(:Lagrangian,Float64,k)
reffe_v = ReferenceFE(:Lagrangian,VectorValue{D,Float64},k)
reffe_λ = ReferenceFE(:Lagrangian,Float64,k)
W = FESpace(model,reffe_w,conformity=:L2)
V = FESpace(model,reffe_v,conformity=:L2)
M = TestFESpace(model_facets,reffe_λ,conformity=:L2,dirichlet_tags="boundary")
Λ = TrialFESpace(M,u)
X = MultiFieldFESpace([V,W])
# Integration machinery
T = Triangulation(model)
∂T = CellBoundary(T)
dT = LebesgueMeasure(T)
d∂T = LebesgueMeasure(∂T)
n = get_normal_vector(∂T)
# FEProblem
a_local((q,u),(w,v)) = ∫( kinv*q⋅v - u*(∇⋅v) - (c*u+q)⋅∇(w) )dT + ∫( ((q+τ*u*n)⋅n)*w )d∂T
l_local_f((w,v)) = ∫( f*w )dT
l_local_m((w,v),m) = ∫( -((c*m -τ*m*n)⋅n)*w - m*(v⋅n) )d∂T
# Variant 1, a_global and l_global in 2 functions
function a_global(λ,μ)
  l_local_η((w,v)) = l_local_m((w,v),λ)
  op_local = LocalHDGOperator(a_local,l_local_η,X,X)
  dq,du = solve(op_local)
  ∫( ((c*η+dq+τ*(du-η)*n)⋅n)*μ )d∂T
end
function l_global(μ)
  op_local = LocalHDGOperator(a_local,l_local_f,X,X)
  q,u = solve(op_local)
  ∫( ((q+τ*u*n)⋅n)*μ )d∂T
end
op = AffineFEOperator(a_global,l_global,Λ,M)
# Variant 2, a_global and l_local in the same function (potential reuse of computation)
function a_l_global(λ,μ)
  l_local_η((w,v)) = l_local_m((w,v),λ)
  op_local = LocalHDGOperator(a_local,l_local_η,l_local_f,X,X)
  (dq,du), (q,u) = solve(op_local)
  ∫( ((c*η+dq+τ*(du-η)*n)⋅n)*μ )d∂T, ∫( ((q+τ*u*n)⋅n)*μ )d∂T
end
op = AffineFEOperator(a_l_global,Λ,M)
# Solve
λ = solve(op)
# Post-process
l_local_η((w,v)) = l_local_m((w,v),λ)
op_local = LocalHDGOperator(a_local,l_local_η,X,X)
q, u = solve(op_local)
# Visualization
E = Triangulation(model_facets)
writevtk(E,"lambda",cellfields=["λ"=>λ])
writevtk(T,"results",cellfields=["q"=>q,"u"=>u])
