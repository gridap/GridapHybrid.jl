
include("test/DarcyRTHTests.jl")

using Profile
using ProfileView
using Gridap
using ExploringGridapHybridization

order=0
domain=(0,1,0,1)
partition=(100,100)
model = CartesianDiscreteModel(domain,partition)
∂T = CellBoundary(model)
DarcyRTHTests.solve_darcy_hybrid_rt(model,∂T,order)
@profile DarcyRTHTests.solve_darcy_hybrid_rt(model,∂T,order)
ProfileView.view()
