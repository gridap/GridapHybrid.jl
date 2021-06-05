module IntegrationTests

using Gridap
using FillArrays
using Gridap.Geometry
using ExploringGridapHybridization

u(x) = VectorValue(2*x[1],x[1]+x[2])
Gridap.divergence(::typeof(u)) = (x) -> 3
p(x) = x[1]-x[2]
∇p(x) = VectorValue(1,-1)
Gridap.∇(::typeof(p)) = ∇p
f(x) = u(x) + ∇p(x)
# Normal component of u(x) on Neumann boundary
function g(x)
  tol=1.0e-14
  if (abs(x[2])<tol)
    return -x[1]-x[2]
  elseif (abs(x[2]-1.0)<tol)
    return x[1]+x[2]
  end
  Gridap.Helpers.@check false
end


function solve_darcy_rt_hdiv()
  domain = (0,1,0,1)
  partition = (1,2)
  order = 0
  model = CartesianDiscreteModel(domain,partition)
  V = FESpace(model,
              ReferenceFE(raviart_thomas,Float64,order),
              conformity=:Hdiv,dirichlet_tags=[5,6])
  Q = FESpace(model,ReferenceFE(lagrangian,Float64,order); conformity=:L2)
  U = TrialFESpace(V,u)
  P = TrialFESpace(Q)
  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])
  trian = Triangulation(model)
  degree = 2
  dΩ = Measure(trian,degree)
  neumanntags = [7,8]
  btrian = BoundaryTriangulation(model,tags=neumanntags)
  degree = 2*(order+1)
  dΓ = Measure(btrian,degree)
  nb = get_normal_vector(btrian)
  a((u, p),(v, q)) = ∫( u⋅v - (∇⋅v)*p + q*(∇⋅u) )*dΩ
  b(( v, q)) = ∫( v⋅f + q*(∇⋅u))*dΩ - ∫((v⋅nb)*p )*dΓ
  op = AffineFEOperator(a,b,X,Y)
  xh = solve(op)
end


# Geometry part
D=2
domain  = (0,1,0,1)
cells   = (1,2)
model   = CartesianDiscreteModel(domain,cells)
model_Γ = BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))

# Functional part
# To investigate what is needed to have an inf-sup stable triplet
# for the RT-H method
order  = 0
reffeᵤ = ReferenceFE(raviart_thomas,Float64,order)
reffeₚ = ReferenceFE(lagrangian,Float64,order)
reffeₗ = ReferenceFE(lagrangian,Float64,order)

# Define test FESpaces
V = TestFESpace(model  , reffeᵤ; conformity=:L2)
Q = TestFESpace(model  , reffeₚ; conformity=:L2)
M = TestFESpace(model_Γ, reffeₗ; conformity=:L2, dirichlet_tags=[7,8])
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

mhscal=get_fe_basis(M)

trian = Triangulation(model)
degree = 2*(order+1)
dΩ = Measure(trian,degree)

neumanntags  = [5,6]
# TO-DO: neumanntrian = Triangulation(model_Γ,tags=neumanntags) this causes
# dcvΓ=∫(mh*g)*dΓ to fail in change_domain ...
neumanntrian = BoundaryTriangulation(model,tags=neumanntags)
degree = 2*(order+1)
dΓ = Measure(neumanntrian,degree)

dcmΩ=∫( vh⋅uh - (∇⋅vh)*ph + qh*(∇⋅uh) )*dΩ
dcvΩ=∫( vh⋅f + qh*(∇⋅u))*dΩ
dcvΓ=∫(mhscal*g)*dΓ

data_mΩ=Gridap.CellData.get_contribution(dcmΩ,dΩ.quad.trian)
data_vΩ=Gridap.CellData.get_contribution(dcvΩ,dΩ.quad.trian)
data_vΓ=Gridap.CellData.get_contribution(dcvΓ,dΓ.quad.trian)

∂T     = CellBoundary(model)
x,w    = quadrature_evaluation_points_and_weights(∂T,2)


#∫( mh*(uh⋅n) )*dK
@time uh_∂T = restrict_to_cell_boundary(∂T,uh)
@time mh_∂T = restrict_to_cell_boundary(∂T,mh)
@time mh_mult_uh_cdot_n=integrate_mh_mult_uh_cdot_n_low_level(∂T,mh_∂T,uh_∂T,x,w)

#∫( (vh⋅n)*lh )*dK
@time vh_∂T = restrict_to_cell_boundary(∂T,vh)
@time lh_∂T = restrict_to_cell_boundary(∂T,lh)
@time vh_cdot_n_mult_lh=integrate_vh_cdot_n_mult_lh_low_level(∂T,vh_∂T,lh_∂T,x,w)

cmat=lazy_map(Broadcasting(+),
              lazy_map(Broadcasting(-),vh_cdot_n_mult_lh,mh_mult_uh_cdot_n),
              data_mΩ)

cvec=data_vΩ

# cell=2
# A11=vcat(hcat(cmat[cell][1,1],cmat[cell][1,2]),hcat(cmat[cell][2,1],0.0))
# A12=vcat(cmat[cell][1,3],zeros(1,4))
# A21=hcat(cmat[cell][3,1],zeros(4))
# S22=-A21*inv(A11)*A12
# b1=vcat(cvec[cell][1],cvec[cell][2])
# y2=-A21*inv(A11)*b1

k=StaticCondensationMap([1,2],[3])
cmat_cvec_condensed=lazy_map(k,cmat,cvec)

fdofsn=lazy_map(Gridap.Arrays.Reindex(get_cell_dof_ids(M)),
                                      get_cell_to_bgcell(dΓ.quad.trian.face_trian))
fdofscb=restrict_facet_dof_ids_to_cell_boundary(∂T,get_cell_dof_ids(M))
assem = SparseMatrixAssembler(M,L)
@time A,b=assemble_matrix_and_vector(assem,(([cmat_cvec_condensed], [fdofscb], [fdofscb]),
                                      ([],[],[]),
                                      ([data_vΓ],[fdofsn])))

x     = A\b
lh    = FEFunction(L,x)

k=BackwardStaticCondensationMap([1,2],[3])
lhₖ= lazy_map(Gridap.Fields.Broadcasting(Gridap.Fields.PosNegReindex(
                      Gridap.FESpaces.get_free_dof_values(lh),lh.dirichlet_values)),
                      fdofscb)
uhphlhₖ=lazy_map(k,cmat,cvec,lhₖ)


lhₑ=lazy_map(Gridap.Fields.BlockMap(3,3),ExploringGridapHybridization.convert_cell_wise_dofs_array_to_facet_dofs_array(∂T,
      lhₖ,get_cell_dof_ids(M)))

assem = SparseMatrixAssembler(Y,X)
lhₑ_dofs=get_cell_dof_ids(X,Triangulation(model_Γ))

uhph_dofs=get_cell_dof_ids(X,Triangulation(model))
uhph_dofs = lazy_map(Gridap.Fields.BlockMap(2,[1,2]),uhph_dofs.args[1],uhph_dofs.args[2])

uh=lazy_map(x->x[1],uhphlhₖ)
ph=lazy_map(x->x[2],uhphlhₖ)
uhphₖ=lazy_map(Gridap.Fields.BlockMap(2,[1,2]),uh,ph)

free_dof_values=assemble_vector(assem,([lhₑ,uhphₖ],[lhₑ_dofs,uhph_dofs]))
xh=FEFunction(X,free_dof_values)

# cell=1
# A11=vcat(hcat(cmat[cell][1,1],cmat[cell][1,2]),hcat(cmat[cell][2,1],0.0))
# A12=vcat(cmat[cell][1,3],zeros(1,4))
# A21=hcat(cmat[cell][3,1],zeros(4))
# b1=vcat(cvec[cell][1],cvec[cell][2])
# x2=lhₖ[cell]
# x1=inv(A11)*(b1-A12*x2)



end # module
