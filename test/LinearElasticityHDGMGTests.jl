# module LinearElasticityHDGMGTests

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
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))
σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε

function A(σ) # compliance tensor
  (1.0 + ν) / E * σ - (ν / E) * tr(σ) * one(σ)
end

# function A(σ) # compliance tensor
#   σ
# end

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

# function invA(ε)
#   ε
# end

function u_exact(x) # Analytical smooth displacement vector
  u1=10.0*sin(pi*x[1])*(1-x[1])*(x[2]-x[2]*x[2])*(1.0-0.5*x[2])
  u2=0.0
  VectorValue(u1,u2)
end

# function u_exact(x) # Linear exact displacement vector
#   VectorValue(x[1] + 2 * x[2], 2 * x[1] + x[2])
# end

function σ_exact(x)
  σ=invA(0.5 * (∇(u_exact)(x) + transpose(∇(u_exact)(x))))
  D=size(σ)[1]
  vals = Gridap.TensorValues._flatten_upper_triangle(reshape(collect(σ.data),(D,D)), Val(D))
  Gridap.TensorValues.SymTensorValue(vals)
end

function f(x)
  (∇ ⋅ (σ_exact))(x)
end

include("P_m.jl")
include("mg_tools.jl")

function build_linear_elasticity_lhdg(model,order)
    # Geometry
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
    # alpha*order*order*(1.0/cells[1]) # HDG stab parameter
    τ = (Float64(num_cells(model)))^(1.0/Float64(D))
    println("τ=$(τ)")

    degree = 2 * (order + 1) # TO-DO: To think which is the minimum degree required
    dΩ = Measure(Ω, degree)
    n = get_cell_normal_vector(∂K)
    d∂K = Measure(∂K, degree)

    # a((σh, uh, uhΓ), (v, ω, μ)) = ∫(v ⊙ (A∘σh) + (∇⋅v) ⋅ uh - ω⋅(∇⋅σh))dΩ -
    #                   ∫((v ⋅ n) ⋅ uhΓ)d∂K +
    #                   #-∫(ω*(σh⋅n-τ*(uh-uhΓ)))*d∂K
    #                   ∫(τ * (ω ⋅ Pₘ(uh, uhΓ, μ, d∂K)))d∂K -
    #                   ∫(τ * (ω ⋅ uhΓ))d∂K +
    #                   #∫(μ*(σh⋅n-τ*(uh-uhΓ)))*d∂K
    #                   ∫(μ ⋅ (σh ⋅ n))d∂K -
    #                   ∫(τ * μ ⋅ Pₘ(uh, uhΓ, μ, d∂K))d∂K +
    #                   ∫(τ * μ ⋅ uhΓ)d∂K
    a((σh, uh, uhΓ), (v, ω, μ)) = ∫(v ⊙ (A∘σh) + (∇⋅v) ⋅ uh - ω⋅(∇⋅σh))dΩ -
                      ∫((v ⋅ n) ⋅ uhΓ)d∂K +
                      #-∫(ω*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                      ∫(τ * (ω ⋅ uh))d∂K -
                      ∫(τ * (ω ⋅ uhΓ))d∂K +
                      #∫(μ*(σh⋅n-τ*(uh-uhΓ)))*d∂K
                      ∫(μ ⋅ (σh ⋅ n))d∂K -
                      ∫(τ * μ ⋅ uh)d∂K +
                      ∫(τ * μ ⋅ uhΓ)d∂K
    l((v, ω, μ)) = ∫(-ω ⋅ f) * dΩ


    op = HybridAffineFEOperator((u, v) -> (a(u, v), l(v)), X, Y, [1, 2], [3])
  end

  function solve_linear_elasticity_lhdg(n,order)
    partition = (0, 1, 0, 1)
    cells = (n, n)
    model = CartesianDiscreteModel(partition, cells)
    D = num_cell_dims(model)
    Ω = Triangulation(ReferenceFE{D}, model)

    op = build_linear_elasticity_lhdg(model, order)

    order_M0 = 1
    reffe_M0 = ReferenceFE(lagrangian, VectorValue{D,Float64}, order_M0; space=:Q)
    M0 = TestFESpace(Ω, reffe_M0; conformity=:H1, dirichlet_tags="boundary")

    degree = 2*order_M0+1
    dΩ = Measure(Ω, degree)
    ∂K = GridapHybrid.Skeleton(model)
    d∂K = Measure(∂K, degree)

    a(u, v) = ∫(ε(v)⊙((A∘σ∘ε(u))))dΩ
    A0 = assemble_matrix(a, M0, M0)

    x1 = zeros(length(op.skeleton_op.op.vector))
    num_iter=nonnested_mg_2_level_v_cycle!(x1, op, A0, M0, dΩ, d∂K;
                                           rtol=1.0e-6, maxiter=1000, smooth_iter=8)

    xexact=op.skeleton_op.op.matrix\op.skeleton_op.op.vector
    num_iter,norm(xexact-x1)
  end

  iters=Int[]
  for n in (5,10,15,20,30,40,50)
     num_iter,_=solve_linear_elasticity_lhdg(n,0)
     append!(iters,num_iter)
  end
  println(iters)


# end # module
