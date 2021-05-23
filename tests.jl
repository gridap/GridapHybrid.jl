using Gridap
using FillArrays
using Gridap.Geometry

include("GridapOverloads.jl")
include("CellBoundary.jl")
include("DensifyInnerMostBlockLevel.jl")

x=rand(3)
y=Array{Vector{Float64},2}(undef,(1,4))
y[1,1]=x
y[1,2]=x
y[1,3]=x
y[1,4]=x
touched  = Array{Bool,2}(undef,(1,4))
touched .= true
cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevel(),Gridap.Fields.ArrayBlock(y,touched))
@test size(cache) == (3,4)
Gridap.Arrays.return_value(cache,DensifyInnerMostBlockLevel(),Gridap.Fields.ArrayBlock(y,touched))

x=rand(1,3)
y=Array{Matrix{Float64}}(undef,(4,))
y[1]=x
y[2]=x
y[3]=x
y[4]=x
touched  = Array{Bool,1}(undef,(4,))
touched .= true
cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevel(),Gridap.Fields.ArrayBlock(y,touched))
Gridap.Arrays.return_value(cache,DensifyInnerMostBlockLevel(),Gridap.Fields.ArrayBlock(y,touched))

x=rand(2,3)
y=Array{Matrix{Float64}}(undef,(3,2))
y[1,1]=x
y[2,1]=x
y[3,1]=x
y[1,2]=x
y[2,2]=x
y[3,2]=x
touched  = Array{Bool,2}(undef,(3,2))
touched .= true
cache=Gridap.Arrays.return_cache(DensifyInnerMostBlockLevel(),Gridap.Fields.ArrayBlock(y,touched))
Gridap.Arrays.return_value(cache,DensifyInnerMostBlockLevel(),Gridap.Fields.ArrayBlock(y,touched))


# Geometry part
D=2
domain  = (0,1,0,1)
cells   = (1,2)
model   = CartesianDiscreteModel(domain,cells)
model_Γ = BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))


# Functional part

# To investigate what is needed to have an inf-sup stable triplet
# for the LDG-H method
order=1
reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
reffeₚ = ReferenceFE(lagrangian,Float64,order)
reffeₗ = ReferenceFE(lagrangian,Float64,order)

# Define test FESpaces
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
_,qh,_ = yh

xh = get_trial_fe_basis(X)
_,_,lh = xh

∂T     = CellBoundary(model)
nowner = get_cell_owner_normal_vector(∂T)
n      = get_cell_normal_vector(∂T)
x,w    = quadrature_evaluation_points_and_weights(∂T,2)

qh_∂T = restrict_to_cell_boundary(∂T,qh)
lh_∂T = restrict_to_cell_boundary(∂T,lh)

i=integrate_low_level(∂T, qh_∂T, lh_∂T, x, w)
