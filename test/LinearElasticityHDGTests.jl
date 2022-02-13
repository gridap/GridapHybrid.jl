module LinearElasticityHDGTests

# Implements HDG formulation proposed in
# https://www.ams.org/journals/mcom/2018-87-309/S0025-5718-2017-03249-X/

using Test
using Gridap
using FillArrays
using Gridap.Geometry
using GridapHybrid
using Plots


include("tmp_fixes_symmetric_valued_lagrangian_reffes.jl")

const ν = 0.3 # Poisson ratio
const E = 1.0 # Young's modulus
function A(σ) # compliance tensor
  (1.0 + ν) / E * σ - (ν / E) * tr(σ) * one(σ)
end
# computes the inverse of the compliance tensor
# applied to the strain tensor
function invA(ε)
  n = prod(size(ε))
  A = zeros(n, n)
  b = zeros(n)
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
  vals = Gridap.TensorValues._flatten_upper_triangle(reshape(A \ b, (D, D)), Val(D))
  Gridap.TensorValues.SymTensorValue(vals)
end

function A(σ) # compliance tensor
  σ
end

function invA(ε) # compliance tensor
  ε
end

function u_exact(x) # Linear exact displacement vector
  VectorValue(x[1] + 2 * x[2], 2 * x[1] + x[2])
end

function u_exact(x) # Analytical smooth displacement vector
  u1=10.0*sin(pi*x[1])*(1-x[1])*(x[2]-x[2]*x[2])*(1.0-0.5*x[2])
  u2=0.0
  VectorValue(u1,u2)
end

function σ_exact(x)
  σ=0.5 * (∇(u_exact)(x) + transpose(∇(u_exact)(x)))
  D=size(σ)[1]
  vals = Gridap.TensorValues._flatten_upper_triangle(reshape(collect(σ.data),(D,D)), Val(D))
  Gridap.TensorValues.SymTensorValue(vals)
end

function f(x)
  (∇ ⋅ (σ_exact))(x)
end

t = Gridap.TensorValues.SymTensorValue{2,Int64,3}(1, 2, 3)
@test t .≈ invA(A(t))
@test t .≈ A(invA(t))

function Gridap.CellData.get_triangulation(f::Gridap.MultiField.MultiFieldCellField)
  s1 = first(f.single_fields)
  trian = get_triangulation(s1)
  trian
end


function solve_linear_elasticity_hdg_symm_tensor(cells,order)
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
    tensor_type = Gridap.Fields.TensorValue{D,D,Float64,D * D}

    # Stress Tensor Space
    reffeσ = ReferenceFE(lagrangian, sym_tensor_type, order; space = :P)
    # Displacement Space
    reffeu = ReferenceFE(lagrangian, VectorValue{D,Float64}, order + 1; space = :P)
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
    τ = 1.0 # HDG stab parameter

    degree = 2 * (order + 1)
    dΩ = Measure(Ω, degree)
    n = get_cell_normal_vector(∂K)
    nₒ = get_cell_owner_normal_vector(∂K)
    d∂K = Measure(∂K, degree)

    yh = get_fe_basis(Y)
    xh = get_trial_fe_basis(X)

    a((σh, uh, uhΓ), (v, ω, μ)) = ∫(v ⊙ (A∘σh) + (∇ ⋅ v) ⋅ uh + ∇(ω) ⊙ σh)dΩ -
                                  ∫((v ⋅ n) ⋅ uhΓ)d∂K -
                                  #-∫(ω*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                                  ∫(ω ⋅ (σh ⋅ n))d∂K +
                                  ∫(τ * (ω ⋅ uh))d∂K -
                                  ∫(τ * (ω ⋅ uhΓ))d∂K +
                                  #∫(μ*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                                  ∫(μ ⋅ (σh ⋅ n))d∂K -
                                  ∫(τ * μ ⋅ uh)d∂K +
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
    #@test sqrt(sum(∫(eu ⋅ eu)dΩ)) < 1.0e-12
    #@test sqrt(sum(∫(eσ ⊙ eσ)dΩ)) < 1.0e-12
  end

  function conv_test(ns,order)
    el2σ = Float64[]
    el2u = Float64[]
    hs = Float64[]
    for n in ns
      l2σ, l2u = solve_linear_elasticity_hdg_symm_tensor((n,n),order)
      h = 1.0/n
      push!(el2σ,l2σ)
      push!(el2u,l2u)
      push!(hs,h)
    end
    (el2σ, el2u, hs)
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x), 1), x) \ y
    linreg[2]
  end

  el2σ, el2u, hs = conv_test([8,16,32,64,128],1)
  plot(hs,[el2σ el2u],
    xaxis=:log, yaxis=:log,
    label=["L2σ" "L2u"],
    shape=:auto,
    xlabel="h",ylabel="error norm")

  println("Slope L2-norm stress: $(slope(hs,el2σ))")
  println("Slope L2-norm      u: $(slope(hs,el2u))")

end # module

# DEBUG USEFUL COMMANDS!

# Γ2=Gridap.Geometry.BoundaryTriangulation(model,[1,2,3,4])
# # n2=get_normal_vector(Γ2)
# dΓ2 = Measure(Γ2, degree)
# dc=∫((v ⋅ n2) ⋅ u_exact)dΓ2

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
