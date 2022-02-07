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

u(x) = VectorValue(1+x[1],1+x[2])
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
W = TestFESpace(Ω  , reffeu; conformity=:L2)
M = TestFESpace(Γ,
                reffe_hat_u;
                conformity=:L2,
                dirichlet_tags=collect(5:8))
Y = MultiFieldFESpace([V,W,M])

# Define trial FEspaces
U = TrialFESpace(V)
P = TrialFESpace(W)
L = TrialFESpace(M,u)
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

σh,uh,uhΓ = xh
v,ω,μ = yh

v⊙(A∘σh)
# (∇⋅v)    # Fails; kills the Julia REPL. To investigate.
# (∇⋅v)⋅uh # Fails; kills the Julia REPL. To investigate.
∇(ω)⊙σh
(v⋅n)⋅uhΓ
ω⋅(σh⋅n)
τ*(ω⋅uh)
# τ*(ω⋅uhΓ) # Fails, issue https://github.com/gridap/GridapHybrid.jl/issues/5
μ⋅(σh⋅n)
#τ*μ⋅uh # Fails, issue https://github.com/gridap/GridapHybrid.jl/issues/5
τ*μ⋅uhΓ
a((σh,uh,uhΓ),(v,ω,μ)) = ∫( v⊙(A∘σh) + (∇⋅v)⋅uh + ∇(ω)⊙σh )dΩ -
                           ∫((v⋅n)⋅uhΓ)d∂K -
                           #-∫(ω*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                           ∫(ω⋅(σh⋅n))d∂K +
                           ∫(τ*(ω⋅uh))d∂K -
                           ∫(τ*(ω⋅uhΓ))d∂K +
                           #∫(μ*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                           ∫(μ⋅(σh⋅n))d∂K -
                           ∫(τ*μ⋅uh)d∂K +
                           ∫(τ*μ⋅uhΓ)d∂K
l((v,ω,μ)) = -∫(ω⋅f)*dΩ

op=HybridAffineFEOperator((u,v)->(a(u,v),l(v)), X, Y, [1,2], [3])
xh=solve(op)

uh,_=xh
e = u -uh
@test sqrt(sum(∫(e⋅e)dΩ)) < 1.0e-12

end # module
