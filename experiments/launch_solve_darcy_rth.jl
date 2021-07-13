
using DrWatson

include("helpers.jl")
include("define_params.jl")
include("driver.jl")

for params in dicts_solve_darcy_rth
  GC.gc()
  run_solve_darcy_rth(params)
end
