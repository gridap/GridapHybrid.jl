module LinearElasticityHDGTests

# Implements HDG formulation proposed in
# https://www.ams.org/journals/mcom/2018-87-309/S0025-5718-2017-03249-X/

using Test
using Gridap
using FillArrays
using Gridap.Geometry
using GridapHybrid
using Plots

const ν = 0.3 # Poisson ratio
const E = 1.0 # Young's modulus
function A(σ) # compliance tensor
  (1.0 + ν) / E * σ - (ν / E) * tr(σ) * one(σ)
end
# computes the inverse of the compliance tensor
# applied to the strain tensor
function invA(ε)
  n = prod(size(ε))
  T=eltype(ε)
  A = zeros(T,n, n)
  b = zeros(T,n)
  for i in eachindex(ε)
    b[i] = ε[i]
  end
  for i = 1:size(A, 1)
    A[i, i] = (1.0 + ν) / E
  end
  D = size(ε)[1]
  for j = 1:D+1:n
    for i = 1:D+1:n
      A[i, j] += -ν / E
    end
  end
  Gridap.TensorValues.TensorValue(reshape(A\b,(D,D)))
end

# function A(σ) # compliance tensor
#   σ
# end

# function invA(ε) # compliance tensor
#   ε
# end

function u_exact(x) # Linear exact displacement vector
  VectorValue(x[1] + 2 * x[2], 2 * x[1] + x[2])
end

function u_exact(x) # Analytical smooth displacement vector
  u1=10.0*sin(pi*x[1])*(1-x[1])*(x[2]-x[2]*x[2])*(1.0-0.5*x[2])
  u2=0.0
  VectorValue(u1,u2)
end

function σ_exact(x)
  σ=invA(0.5 * (∇(u_exact)(x) + transpose(∇(u_exact)(x))))
  D=size(σ)[1]
  vals = Gridap.TensorValues._flatten_upper_triangle(reshape(collect(σ.data),(D,D)), Val(D))
  Gridap.TensorValues.SymTensorValue(vals)
end

function f(x)
  (∇ ⋅ (σ_exact))(x)
end

t = Gridap.TensorValues.SymTensorValue{2,Float64,3}(1, 2, 3)
@test t .≈ invA(A(t))
@test t .≈ A(invA(t))

include("P_m.jl")

function solve_linear_elasticity_hdg_symm_tensor(
    cells,order;alpha=1.0,bulk_to_skeleton_projection=true,write_results=false)
    # Geometry
    partition = (0,1,0,1)
    model = CartesianDiscreteModel(partition, cells)
    D = num_cell_dims(model)
    Ω = Triangulation(ReferenceFE{D}, model)
    Γ = Triangulation(ReferenceFE{D-1}, model)
    ∂K = GridapHybrid.Skeleton(model)

    # Reference FEs
    num_comp = D * (D + 1) ÷ 2 # Number of components of a symmetric tensor in D-dim space
    sym_tensor_type = Gridap.Fields.SymTensorValue{D,Float64,num_comp}

    # Stress Tensor Space
    reffeσ = ReferenceFE(lagrangian, sym_tensor_type, order; space = :P)
    # Displacement Space
    reffeu = ReferenceFE(lagrangian, VectorValue{D,Float64}, order+1; space = :P)
    # Trace of Displacements Space
    reffe_hat_u = ReferenceFE(lagrangian, VectorValue{D,Float64}, order; space = :P)

    # Define test FESpaces
    V = TestFESpace(Ω, reffeσ; conformity = :L2)
    W = TestFESpace(Ω, reffeu; conformity = :L2)
    M = TestFESpace(Γ,
                    reffe_hat_u;
                    conformity = :L2,
                    dirichlet_tags = collect(5:8))
    Y = MultiFieldFESpace([V, W, M])

    # Define trial FEspaces
    U = TrialFESpace(V)
    P = TrialFESpace(W)
    L = TrialFESpace(M, u_exact)
    X = MultiFieldFESpace([U, P, L])

    # FE formulation params
    τ = cells[1] # alpha*order*order*(1.0/cells[1]) # HDG stab parameter
    println("h=$(1.0/cells[1]) τ=$(τ)")

    degree = 2 * (order + 1) # TO-DO: To think which is the minimum degree required
    dΩ = Measure(Ω, degree)
    n = get_cell_normal_vector(∂K)
    nₒ = get_cell_owner_normal_vector(∂K)
    d∂K = Measure(∂K, degree)

    yh = get_fe_basis(Y)
    xh = get_trial_fe_basis(X)

    function project(uh, uhΓ, μ, d∂K; bulk_to_skeleton_projection=true)
       if (bulk_to_skeleton_projection)
         Pₘ(uh, uhΓ, μ, d∂K)
       else
         uh
       end
    end

    # Testing for correctness of Pₘ projection
    (v, ω, μ) = yh
    (σh, uh, uhΓ) = xh
    dc=∫(μ⋅(uh-Pₘ(uh, uhΓ, μ, d∂K)))d∂K
    for k in dc.dict[d∂K.quad.trian]
       @test all(k[3,2] .< 1.0e-14)
    end

    a((σh, uh, uhΓ), (v, ω, μ)) = ∫(v ⊙ (A∘σh) + (∇⋅v) ⋅ uh - ω⋅(∇⋅σh))dΩ -
                      ∫((v ⋅ n) ⋅ uhΓ)d∂K +
                      #-∫(ω*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                      ∫(τ * (ω ⋅ project(uh, uhΓ, μ, d∂K;bulk_to_skeleton_projection=bulk_to_skeleton_projection)))d∂K -
                      ∫(τ * (ω ⋅ uhΓ))d∂K +
                      #∫(μ*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                      ∫(μ ⋅ (σh ⋅ n))d∂K -
                      ∫(τ * μ ⋅ project(uh, uhΓ, μ, d∂K;
                      bulk_to_skeleton_projection=bulk_to_skeleton_projection))d∂K +
                      ∫(τ * μ ⋅ uhΓ)d∂K
    l((v, ω, μ)) = ∫(-ω ⋅ f) * dΩ


    op = HybridAffineFEOperator((u, v) -> (a(u, v), l(v)), X, Y, [1, 2], [3])
    xh = solve(op)

    σh, uh, uhΓ = xh
    eσ = σ_exact - σh
    eu = u_exact - uh
    norm2_σ=sqrt(sum(∫(eσ ⊙ eσ)dΩ))
    norm2_u=sqrt(sum(∫(eu ⋅ eu)dΩ))
    norm2_σ,norm2_u

    if (write_results)
      writevtk(Ω,"results_$(cells)_k=$(order)",cellfields=["σh"=>σh,"uh"=>uh,"σ_exact"=>σ_exact, "error"=>eσ])
    end
    return     norm2_σ,norm2_u
    #@test sqrt(sum(∫(eu ⋅ eu)dΩ)) < 1.0e-12
    #@test sqrt(sum(∫(eσ ⊙ eσ)dΩ)) < 1.0e-12
  end


  function conv_test(ns,order;alpha=1.0,bulk_to_skeleton_projection=false)
    el2σ = Float64[]
    el2u = Float64[]
    hs = Float64[]
    for n in ns
      l2σ, l2u = solve_linear_elasticity_hdg_symm_tensor((n,n),order;
                                                          alpha=alpha,bulk_to_skeleton_projection=bulk_to_skeleton_projection)
      println(l2σ, " ", l2u)
      h = 1.0/n
      push!(el2σ,l2σ)
      push!(el2u,l2u)
      push!(hs,h)
    end
    println(el2σ)
    println(el2u)
    (el2σ, el2u, hs)
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x), 1), x) \ y
    linreg[2]
  end

  ns=[8,16,32,64,128]
  order=1
  el2σ_noPM, el2u_noPM, hs = conv_test(ns,order;bulk_to_skeleton_projection=false)
  slopek  =[Float64(ni)^(-(order)) for ni in ns]
  slopekp1=[Float64(ni)^(-(order+1)) for ni in ns]
  slopekp2=[Float64(ni)^(-(order+2)) for ni in ns]
  plot(hs,[el2σ_noPM el2u_noPM slopek slopekp1 slopekp2],
    xaxis=:log, yaxis=:log,
    label=["L2σ (measured)" "L2u (measured)" "slope k" "slope k+1" "slope k+2"],
    shape=:auto,
    xlabel="h",ylabel="L2 error (No PM)",legend=:bottomright)

  println("Slope L2-norm stress (no PM): $(slope(hs,el2σ_noPM))")
  println("Slope L2-norm      u (no PM): $(slope(hs,el2u_noPM))")

  el2σ_PM, el2u_PM, hs = conv_test([8,16,32,64,128],1;alpha=1.0,bulk_to_skeleton_projection=true)
  plot(hs,[el2σ_PM el2u_PM slopek slopekp1 slopekp2],
    xaxis=:log, yaxis=:log,
    label=["L2σ (measured)" "L2u (measured)" "slope k" "slope k+1" "slope k+2"],
    shape=:auto,
    xlabel="h",ylabel="L2 error norm (PM)",legend=:bottomright)

  println("Slope L2-norm stress (PM): $(slope(hs,el2σ_PM))")
  println("Slope L2-norm      u (PM): $(slope(hs,el2u_PM))")

end # module

# # # Geometry
# partition = (0,2,0,1)
# cells=(2,1)
# model = CartesianDiscreteModel(partition, cells)
# D = num_cell_dims(model)
# Ω = Triangulation(ReferenceFE{D}, model)
# Γ = Triangulation(ReferenceFE{D-1}, model)
# ∂K = GridapHybrid.Skeleton(model)

# # Stress Tensor Space
# order=1
# num_comp = D * (D + 1) ÷ 2 # Number of components of a symmetric tensor in D-dim space
# sym_tensor_type = Gridap.Fields.SymTensorValue{D,Float64,num_comp}
# reffeσ = ReferenceFE(lagrangian, sym_tensor_type, order; space = :P)
# # Displacement Space
# reffeu = ReferenceFE(lagrangian, VectorValue{D,Float64}, order+1; space = :P)
# # Trace of Displacements Space
# reffe_hat_u = ReferenceFE(lagrangian, VectorValue{D,Float64}, order; space = :P)

# # Define test FESpaces
# V = TestFESpace(Ω, reffeσ; conformity = :L2)
# W = TestFESpace(Ω, reffeu; conformity = :L2)
# M = TestFESpace(Γ,
#                 reffe_hat_u;
#                 conformity = :L2,
#                 dirichlet_tags = collect(5:8))
# Y = MultiFieldFESpace([V, W, M])

# # Define trial FEspaces
# U = TrialFESpace(V)
# P = TrialFESpace(W)
# L = TrialFESpace(M, u_exact)
# X = MultiFieldFESpace([U, P, L])


# yh = get_fe_basis(Y)
# xh = get_trial_fe_basis(X)


# degree = 2 * (order + 1)
# dΩ = Measure(Ω, degree)
# n = get_cell_normal_vector(∂K)
# nₒ = get_cell_owner_normal_vector(∂K)
# d∂K = Measure(∂K, degree)

# v, ω, μ = yh
# σh, uh, uhΓ = xh
# Pmuh=Pₘ(uh, uhΓ, μ, d∂K)
# x=get_cell_points(d∂K.quad)
# pmuh_at_x=Pmuh(x)
# # DEBUG USEFUL COMMANDS!

# Γ2=Gridap.Geometry.BoundaryTriangulation(model,collect(1:7))
# # n2=get_normal_vector(Γ2)
# dΓ2 = Measure(Γ2, degree)
# #dc=∫((v ⋅ n2) ⋅ u_exact)dΓ2

# M2 = TestFESpace(Γ2,
#                  reffe_hat_u;
#                  conformity = :L2)
# L2 = TrialFESpace(M2)

# uh2=get_fe_basis(W)
# function f(c::Integer, i::Integer,x::Point)
#   lazy_map(evaluate,Gridap.CellData.get_data(uh2),[[x],[x]])[c][1,i]
# end
# n=length(Gridap.CellData.get_data(uh2)[1])
# x2=get_cell_points(dΓ2.quad)
# cfids=[[1,2,3,4],[5,6,4,7]]

# for c=1:2
#   for i=1:n
#     function f(x)
#      f(c,i,x)
#     end
#     a(u,v)=∫(v⋅u)dΓ2
#     l(v)=∫(v⋅f)dΓ2
#     op=AffineFEOperator(a,l,M2,L2)
#     phi_i_projected=solve(op)
#     a=phi_i_projected(x2)

#     println(a)

#     for f=1:4
#       println("$(c) $(i) $(f)")
#       r=a[cfids[c][f]]
#       s=pmuh_at_x[c][f][1,2][1][:,1,i]
#       println(r)
#       println(s)
#       @assert all(r .≈ s)
#       println("")
#     end
#   end
# end


# M2 = TestFESpace(Γ2,
#   reffe_hat_u;
#   conformity = :L2,
#   dirichlet_tags = collect(5:8))
# μ2 = get_fe_basis(M2)
# dc=∫(τ * μ2 ⋅ uh_exact)dΓ2
# dc2=∫(τ * μ ⋅ uh_exact)d∂K


# uh_exact=interpolate_everywhere(u_exact,W)


# r((v, ω, μ)) = ∫(μ ⋅ (σ_exact ⋅ n))d∂K - ∫(τ * μ ⋅ uh_exact)d∂K + ∫(τ * μ ⋅ uh_exact)d∂K
# r((v, ω, μ)) = ∫(τ * μ ⋅ uh_exact)d∂K
# rvec=assemble_vector(r,Y)

# dc=∫(μ ⋅ (σ_exact ⋅ n))d∂K - ∫(τ * μ ⋅ uh_exact)d∂K + ∫(τ * μ ⋅ uh_exact)d∂K

# dc.dict[keys(dc.dict)...]

# dc=∫(v ⊙ (A∘σ_exact) + (∇⋅v) ⋅ uh_exact)dΩ - ∫((v ⋅ n) ⋅ uh_exact)d∂K
# dc.dict[x[1]][1][1]
# dc.dict[x[2]][1][1]

# dc.dict[x[1]][1][1]+dc.dict[x[2]][1][1]

# r((v, ω, μ)) = ∫(∇(ω) ⊙ σ_exact)dΩ -
#                ∫(ω ⋅ (σ_exact ⋅ n))d∂K + ∫(τ * (ω ⋅ uh_exact))d∂K - ∫(τ * (ω ⋅ uh_exact))d∂K

# dc=∫(∇(ω) ⊙ σ_exact)dΩ -
#    ∫(ω ⋅ (σ_exact ⋅ n))d∂K + ∫(τ * (ω ⋅ uh_exact))d∂K - ∫(τ * (ω ⋅ uh_exact))d∂K

# rvec=assemble_vector(r,Y)

# x=[]
# for k in keys(dc.dict)
#   push!(x,k)
# end

# dc.dict[x[2]][1][2]+dc.dict[x[1]][1][2]

# r((v, ω, μ)) = ∫(v ⊙ (A∘σ_exact) + (∇⋅v) ⋅ u_exact)dΩ - ∫((v ⋅ n) ⋅ u_exact)d∂K
# r1((v, ω, μ)) = ∫(v ⊙ (A∘σ_exact))dΩ # ∫((v ⋅ n) ⋅ u_exact)d∂K
# r2((v, ω, μ)) = ∫((∇⋅v) ⋅ u_exact)dΩ # ∫((v ⋅ n) ⋅ u_exact)d∂K
# r3((v, ω, μ)) = ∫((v ⋅ n) ⋅ u_exact)d∂K
# r1vec=assemble_vector(r,Y)
# r2vec=assemble_vector(r2,Y)
# r3vec=assemble_vector(r3,Y)

# tensor_type = Gridap.Fields.TensorValue{D,D,Float64,D * D}
# # # Stress Tensor Space
# reffeσ1 = ReferenceFE(lagrangian, Float64, order; space = :P)

# # Define test FESpaces
# V1 = TestFESpace(Ω, reffeσ1; conformity = :L2)

# v1 = get_fe_basis(V1)
# ∇⋅v1

# ∇⋅v1

# x = get_cell_points(dΩ.quad)
# div1=(∇⋅v1)(x)

# v = get_fe_basis(V)
# ∇⋅v
# x = get_cell_points(dΩ.quad)
# div=(∇⋅v)(x)


# function vi_at_x(v,i,x)
#   # V is a cell Field
#   # i is the identifier of the basis function
#   lazy_map(evaluate,Gridap.CellData.get_data(v),[x])[1][1][i]
# end

# i=1
# function vi(x)
#   vi_at_x(v,i,x)
# end

# dvi=divergence(vi)
