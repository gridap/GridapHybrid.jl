using DrWatson

jobname(args...) = replace(savename(args...),"="=>"")

function collect_from_dicts(experiment,dicts)
  results = []
  for params in dicts
    outfile = datadir("experiments",jobname(experiment,params,"bson"))
    if isfile(outfile)
      out = load(outfile)
      push!(results,out)
    end
  end
  results
end

function replace_strings_by_symbols(r)
  d = Dict{Symbol,Any}()
  for (k,v) in r
    d[Symbol(k)] = v
  end
  d
end

function compute_min_time(results,ks)
  dict = Dict()
  for out in results
    key = Tuple([ out[k] for k in ks])
    if haskey(dict,key)
      out2 = dict[key]
      for (k,v) in out
        if occursin("[s]",k)
          out2[k] = min(out2[k],v)
        end
      end
    else
      dict[key] = out
    end
  end
  collect(values(dict))
end

function restrict_to_keys(i,ks)
  o = typeof(i)()
  for k in ks
    o[k] = i[k]
  end
  o
end

function run_solve_darcy_rth(params)
  outfile = datadir("experiments",jobname("solve_darcy_rth",params,"bson"))
  if isfile(outfile)
    println("$outfile (done already)")
    return nothing
  end
  print("$outfile (running)")
  n = params[:n]
  k = params[:k]
  d = params[:d]
  dict = Experiments.solve_darcy_rth(n,k,d)
  save(outfile,dict)
  println(" (done)")
  nothing
end

function run_solve_darcy_rt(params)
  outfile = datadir("experiments",jobname("solve_darcy_rt",params,"bson"))
  if isfile(outfile)
    println("$outfile (done already)")
    return nothing
  end
  print("$outfile (running)")
  n = params[:n]
  k = params[:k]
  d = params[:d]
  dict = Experiments.solve_darcy_rt(n,k,d)
  save(outfile,dict)
  println(" (done)")
  nothing
end
