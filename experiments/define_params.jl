using DrWatson

nruns = 10
params_solve_darcy_rth = []

push!(params_solve_darcy_rth,Dict(
 :n => [10,20,30,40,50,60,70,80,90,100],
 #:n => [5,10,15,20,25,30],
 :k => [0,],
 :d => [2,],
 :run => collect(1:nruns)
))

# push!(params_solve_darcy_rth,Dict(
#  :n => [10,20,30,40,50],
#  :k => [2,],
#  :run => collect(1:nruns)
# ))

# push!(params_solve_darcy_rth,Dict(
#  :n => [1,2,3,4,5,6],
#  :k => [1,],
#  :run => collect(1:nruns)
# ))

# push!(params_solve_darcy_rth,Dict(
#  :n => [1,2,3,4,5],
#  :k => [2,],
#  :run => collect(1:nruns)
# ))

dicts_solve_darcy_rth = vcat(map(dict_list,params_solve_darcy_rth)...)
