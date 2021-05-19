using Gridap
using FillArrays
using Gridap.Geometry

include("GridapOverloads.jl")
include("CellBoundary.jl")


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
∂Tq    = quadrature_evaluation_points_and_weights(∂T,2)

qh_cb = restrict_to_cell_boundary(∂T,qh)
lh_cb = restrict_to_cell_boundary(∂T,lh)

qh_cb_x = lazy_map(evaluate,qh_cb,∂Tq[1])
lh_cb_x = lazy_map(evaluate,lh_cb,∂Tq[1])

qh_mult_lh_cb_x = lazy_map(Gridap.Fields.BroadcastingFieldOpMap(*),qh_cb_x,lh_cb_x)
