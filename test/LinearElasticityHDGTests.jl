module LinearElasticityHDGTests

# Implements HDG formulation proposed in
# https://www.ams.org/journals/mcom/2018-87-309/S0025-5718-2017-03249-X/

using Test
using Gridap
using FillArrays
using Gridap.Geometry
using GridapHybrid

include("tmp_fixes_symmetric_valued_lagrangian_reffes.jl")

const ν=0.3 # Poisson ratio
const E=1.0 # Young's modulus
function A(σ) # compliance tensor
  (1.0+ν)/E*σ-(ν/E)*tr(σ)*one(σ)
end

# computes the inverse of the compliance tensor
# applied to the strain tensor
function invA(ε)
  n=prod(size(ε))
  A=zeros(n,n)
  b=zeros(n)
  for i in eachindex(ε)
    b[i]=ε[i]
  end
  for i=1:size(A,1)
    A[i,i]=(1.0+ν)/E
  end
  D=size(ε)[1]
  for j=1:D+1:n
    for i=1:D+1:n
      A[i,j]+=-ν/E
    end
  end
  vals=Gridap.TensorValues._flatten_upper_triangle(reshape(A\b,(D,D)),Val(D))
  Gridap.TensorValues.SymTensorValue(vals)
end

t=Gridap.TensorValues.SymTensorValue{2, Int64, 3}(1, 2, 3)
@test t .≈ invA(A(t))
@test t .≈ A(invA(t))

# u(x) = VectorValue(1+x[1],1+x[2])
# Gridap.divergence(::typeof(u)) = (x) -> 2
# p(x) = -3.14
# ∇p(x) = VectorValue(0,0)
# Gridap.∇(::typeof(p)) = ∇p
# f(x) = u(x) + ∇p(x)
# # Normal component of u(x) on Neumann boundary
# function g(x)
#   tol=1.0e-14
#   if (abs(x[2])<tol)
#     return -x[2] #-x[1]-x[2]
#   elseif (abs(x[2]-1.0)<tol)
#     return x[2] # x[1]+x[2]
#   end
#   Gridap.Helpers.@check false
# end

function Gridap.CellData.get_triangulation(f::Gridap.MultiField.MultiFieldCellField)
  s1 = first(f.single_fields)
  trian = get_triangulation(s1)
  trian
end

partition = (0,1,0,1)
cells = (2,2)
model = CartesianDiscreteModel(partition,cells)

# Geometry
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D},model)
Γ = Triangulation(ReferenceFE{D-1},model)
∂K = GridapHybrid.Skeleton(model)

# Reference FEs
order           = 1
num_comp        = D*(D+1)÷2 # Number of components of a symmetric tensor in D-dim space
sym_tensor_type = Gridap.Fields.SymTensorValue{D,Float64,num_comp}

# Stress Tensor Space
reffeσ      = ReferenceFE(lagrangian,sym_tensor_type,order;space=:P)
# Displacement Space
reffeu      = ReferenceFE(lagrangian,VectorValue{D,Float64},order+1;space=:P)
# Trace of Displacements Space
reffe_hat_u = ReferenceFE(lagrangian,VectorValue{D,Float64},order;space=:P)

# Define test FESpaces
V = TestFESpace(Ω  , reffeσ; conformity=:L2)
Q = TestFESpace(Ω  , reffeu; conformity=:L2)
M = TestFESpace(Γ,
                reffe_hat_u;
                conformity=:L2,
                dirichlet_tags=collect(5:8))
Y = MultiFieldFESpace([V,Q,M])

# Define trial FEspaces
U = TrialFESpace(V)
P = TrialFESpace(Q)
L = TrialFESpace(M,p)
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

(uh,ph,lh) = xh
(vh,qh,mh) = yh

vh_∂K = Gridap.CellData.change_domain(vh,∂K,ReferenceDomain())
mh_∂K = Gridap.CellData.change_domain(mh,∂K,ReferenceDomain())

# Left-to-right evaluation of ∫(τ*(mh*(ph*(n⋅nₒ))))d∂K
# TO DO: currently fails in op1*ph !!!
# op1 = (n⋅nₒ)
# op1 = τ*mh
# op2 = op1*ph
# op3 = op2*op1

# Right-to-left evaluation
# Works
op1=ph*(n⋅nₒ)
mh*op1
τ*(mh*(ph*(n⋅nₒ)))
lh*(n⋅nₒ)

# Left-to-right evaluation
# Works
τ*mh*lh*(n⋅nₒ)

# Right-to-left evaluation
# Works
τ*(mh*(lh*(n⋅nₒ)))

# Evaluation of ∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
# qh*(uh⋅n)
# Works
qh*(uh⋅n)
(qh*uh)⋅n

# τ*(qh*(ph*(n⋅nₒ)))
# Right-to-left
# Works
op1=ph*(n⋅nₒ)
op2=qh*op1
op3=τ*op2

# Left-to-right
# τ*qh*ph*(n⋅nₒ)
# Fails! in op3 TO-DO
op1=τ*qh
op2=op1*ph
#op3=op2*(n⋅nₒ)

# τ*(qh*(lh*(n⋅nₒ)))
# Right-to-left
# Works
op1=lh*(n⋅nₒ)
op2=qh*op1
op3=τ*op2

# Left-to-right
# Fails! in op2 TO-DO
op1=τ*qh
#op2=op1*lh
#op3=op2*(n⋅nₒ)

∫(mh*(uh⋅n))d∂K

a((uh,ph,lh),(vh,qh,mh)) = ∫( vh⋅uh - (∇⋅vh)*ph - ∇(qh)⋅uh )dΩ +
                           ∫((vh⋅n)*lh)d∂K +
                           #∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                           ∫(qh*(uh⋅n))d∂K +
                           ∫(τ*(qh*(ph*(n⋅nₒ))))d∂K - # Force right-to-left evaluation
                           ∫(τ*(qh*(lh*(n⋅nₒ))))d∂K + # Force right-to-left evaluation
                           #∫(mh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                           ∫(mh*(uh⋅n))d∂K +
                           ∫(τ*(mh*(ph*(n⋅nₒ))))d∂K - # Force right-to-left evaluation
                           ∫(τ*mh*lh*(n⋅nₒ))d∂K
l((vh,qh,mh)) = ∫( vh⋅f + qh*(∇⋅u))*dΩ

op=HybridAffineFEOperator((u,v)->(a(u,v),l(v)), X, Y, [1,2], [3])
xh=solve(op)

uh,_=xh
e = u -uh
@test sqrt(sum(∫(e⋅e)dΩ)) < 1.0e-12

end # module
