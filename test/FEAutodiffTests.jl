module FEAutodiffTests

using Gridap
using GridapHybrid
using Test

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

# Geometry
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D},model)
Γ = Triangulation(ReferenceFE{D-1},model)
∂K = GridapHybrid.Skeleton(model)

# Reference FEs
order  = 1
reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order;space=:P)
reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
reffeₗ = ReferenceFE(lagrangian,Float64,order;space=:P)

# Define test FESpaces
V = TestFESpace(Ω  , reffeᵤ; conformity=:L2)
Q = TestFESpace(Ω  , reffeₚ; conformity=:L2)
M = TestFESpace(Γ,reffeₗ;conformity=:L2)
Y = MultiFieldFESpace([V,Q,M])

# Define trial FEspaces
U = TrialFESpace(V)
P = TrialFESpace(Q)
L = TrialFESpace(M)
X = MultiFieldFESpace([U, P, L])

# FE formulation params
τ = 1.0 # HDG stab parameter

degree = 2*(order+1)
dΩ     = Measure(Ω,degree)
n      = get_cell_normal_vector(∂K)
nₒ     = get_cell_owner_normal_vector(∂K)
d∂K    = Measure(∂K,degree)

yh = get_fe_basis(Y)
xh = get_trial_fe_basis(X)
fh = FEFunction(X,ones(num_free_dofs(X)))

res((uh,ph,lh),(vh,qh,mh)) =
                          ∫( vh⋅uh - (∇⋅vh)*ph - ∇(qh)⋅uh )dΩ +
                          ∫((vh⋅n)*lh)d∂K +
                          #∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                          ∫(qh*(uh⋅n))d∂K +
                          ∫(τ*qh*ph*(n⋅nₒ))d∂K -
                          ∫(τ*qh*lh*(n⋅nₒ))d∂K +
                          #∫(mh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                          ∫(mh*(uh⋅n))d∂K +
                          ∫(τ*mh*ph*(n⋅nₒ))d∂K -
                          ∫(τ*mh*lh*(n⋅nₒ))d∂K

jac((uh,ph,lh),(duh,dph,dlh),(vh,qh,mh)) =
                          ∫( vh⋅duh - (∇⋅vh)*dph - ∇(qh)⋅duh )dΩ +
                          ∫((vh⋅n)*dlh)d∂K +
                          #∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                          ∫(qh*(duh⋅n))d∂K +
                          ∫(τ*qh*dph*(n⋅nₒ))d∂K -
                          ∫(τ*qh*dlh*(n⋅nₒ))d∂K +
                          #∫(mh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                          ∫(mh*(duh⋅n))d∂K +
                          ∫(τ*mh*dph*(n⋅nₒ))d∂K -
                          ∫(τ*mh*dlh*(n⋅nₒ))d∂K

dc=GridapHybrid._merge_bulk_and_skeleton_matrix_contributions(jac(fh,xh,yh))
A=assemble_matrix(dc,X,Y)
cell_j = get_array(dc)
cell_j_auto = get_array(jacobian(x->res(x,yh),fh))

function compare(x,y)
  res = size(x) == size(y)
  for i in eachindex(x)
    if x.touched[i] && y.touched[i]
      res = res && (x.array[i]≈y.array[i])
    end
    res
  end
end

@test Gridap.Arrays.test_array(cell_j_auto,cell_j,≈)

end # module
