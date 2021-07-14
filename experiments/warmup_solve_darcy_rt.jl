using Gridap
using DrWatson

include("helpers.jl")
include("define_params.jl")
include("driver.jl")
Experiments.solve_darcy_rt(10,0,2)
Experiments.solve_darcy_rt(10,0,2)
Experiments.solve_darcy_rt(5,0,3)
Experiments.solve_darcy_rt(5,0,3)
