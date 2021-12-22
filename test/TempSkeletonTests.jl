using Gridap
using ExploringGridapHybridization

function Gridap.CellData.get_triangulation(f::Gridap.MultiField.MultiFieldCellField)
  s1 = first(f.single_fields)
  trian = get_triangulation(s1)
  trian
end

partition = (0,1,0,1)
cells = (5,5)
model = CartesianDiscreteModel(partition,cells)


D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D},model)
Γ = Triangulation(ReferenceFE{D-1},model)
∂K = TempSkeleton(model)
println(get_cell_coordinates(∂K))

#v = CellField(x->x[1]+x[2],Ω)
#q = CellField(x->x[1]-x[2]+1,Γ)
#writevtk(Ω,"draft_Ω",cellfields=["v"=>v])
#writevtk(Γ,"draft_Γ",cellfields=["q"=>q])
#writevtk(∂K,"draft_∂K",offset=0.15,cellfields=["q"=>q,"v"=>v])

# Functional part
# To investigate what is needed to have an inf-sup stable triplet
# for the RT-H method
order=0
reffeᵤ = ReferenceFE(raviart_thomas,Float64,order)
reffeₚ = ReferenceFE(lagrangian,Float64,order)
reffeₗ = ReferenceFE(lagrangian,Float64,order)

# Define test FESpaces
V = TestFESpace(Ω  , reffeᵤ; conformity=:L2)
Q = TestFESpace(Ω  , reffeₚ; conformity=:L2)
M = TestFESpace(Γ,
                reffeₗ;
                conformity=:L2,
                dirichlet_tags=collect(5:8))
Y = MultiFieldFESpace([V,Q,M])

# Create trial spaces
U = TrialFESpace(V)
P = TrialFESpace(Q)
L = TrialFESpace(rand(num_free_dofs(M)),M)
X = MultiFieldFESpace([U, P, L])

yh = get_fe_basis(Y)
vh,qh,mh = yh

xh = get_trial_fe_basis(X)
uh,ph,lh = xh

degree=2

n   = get_cell_normal_vector(∂K)
d∂K = Measure(∂K,degree)
(vh⋅n)
dc=∫((vh⋅n)*lh)d∂K
dc_array=Gridap.CellData.get_contribution(dc,∂K)
