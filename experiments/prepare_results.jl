using DrWatson
using DataFrames

include("helpers.jl")
include("define_params.jl")

results = collect_from_dicts("solve_darcy_rth",dicts_solve_darcy_rth)
ks = collect(keys(first(dicts_solve_darcy_rth)))
uks = map(string,ks[findall(k->k!=:run,ks)])
results = compute_min_time(results,uks)
results = map(replace_strings_by_symbols,results)
df_solve_darcy_rth = DataFrame(results)
sort!(df_solve_darcy_rth,[:d,:k,:n])

results = collect_from_dicts("solve_darcy_rt",dicts_solve_darcy_rth)
ks = collect(keys(first(dicts_solve_darcy_rth)))
uks = map(string,ks[findall(k->k!=:run,ks)])
results = compute_min_time(results,uks)
results = map(replace_strings_by_symbols,results)
df_solve_darcy_rt = DataFrame(results)
sort!(df_solve_darcy_rt,[:d,:k,:n])
