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
uh_∂T = restrict_to_cell_boundary(∂T,uh)
mh_∂T = restrict_to_cell_boundary(∂T,mh)
mh_mult_uh_cdot_n=integrate_mh_mult_uh_cdot_n_low_level(∂T,mh_∂T,uh_∂T,x,w)

#∫( (vh⋅n)*lh )*dK
vh_∂T = restrict_to_cell_boundary(∂T,vh)
lh_∂T = restrict_to_cell_boundary(∂T,lh)
vh_cdot_n_mult_lh=integrate_vh_cdot_n_mult_lh_low_level(∂T,vh_∂T,lh_∂T,x,w)

cmat=lazy_map(Broadcasting(+),
              lazy_map(Broadcasting(-),data_mΩ,vh_cdot_n_mult_lh),
              mh_mult_uh_cdot_n)

cvec=data_vΩ

k=StaticCondensationMap([1,2],[3])
cmat_cvec_condensed=lazy_map(k,cmat,cvec)
cmat_condensed=lazy_map(x->x[1],cmat_cvec_condensed)
cvec_condensed=lazy_map(x->x[2],cmat_cvec_condensed)


# Include Neumann term
w = []
r = []
push!(w,data_vΓ)
push!(r,lazy_map(Gridap.Arrays.Reindex(get_cell_dof_ids(M)),
                 get_cell_to_bgcell(dΓ.quad.trian.face_trian)))

# Include cell term
fdofscb=restrict_facet_dof_ids_to_cell_boundary(∂T,get_cell_dof_ids(M))
push!(w,cvec_condensed)
push!(r,fdofscb)

assem = SparseMatrixAssembler(M,L)
b     = assemble_vector(assem,(w,r))
A     = assemble_matrix(assem,( [cmat_condensed], [fdofscb], [fdofscb]))
x     = A\b
xh    = FEFunction(L,x)

k=BackwardStaticCondensationMap([1,2],[3])
xhₖ= lazy_map(Gridap.Fields.Broadcasting(Gridap.Fields.PosNegReindex(
                      Gridap.FESpaces.get_free_dof_values(xh),xh.dirichlet_values)),
                      fdofscb)
x1x2=lazy_map(k,cmat,cvec,xhₖ)
end
