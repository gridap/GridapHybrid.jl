
module DarcyHDGMGTests

  using Gridap
  using GridapHybrid
  using FillArrays
  using IterativeSolvers
  using Printf

  include("P_m.jl")
  include("mg_tools.jl")

  u(x) = VectorValue(1 + x[1], 1 + x[2])
  Gridap.divergence(::typeof(u)) = (x) -> 2
  p(x) = -3.14
  ∇p(x) = VectorValue(0, 0)
  Gridap.∇(::typeof(p)) = ∇p
  f(x) = u(x) + ∇p(x)
  # Normal component of u(x) on Neumann boundary
  function g(x)
    tol = 1.0e-14
    if (abs(x[2]) < tol)
      return -x[2] #-x[1]-x[2]
    elseif (abs(x[2] - 1.0) < tol)
      return x[2] # x[1]+x[2]
    end
    Gridap.Helpers.@check false
  end

  function build_darcy_lhdg(model, order)
    # Geometry
    D = num_cell_dims(model)
    Ω = Triangulation(ReferenceFE{D}, model)
    Γ = Triangulation(ReferenceFE{D - 1}, model)
    ∂K = GridapHybrid.Skeleton(model)

    # Reference FEs
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order; space=:P)
    reffeₚ = ReferenceFE(lagrangian, Float64, order; space=:P)
    reffeₗ = ReferenceFE(lagrangian, Float64, order; space=:P)

    # Define test FESpaces
    V = TestFESpace(Ω, reffeᵤ; conformity=:L2)
    Q = TestFESpace(Ω, reffeₚ; conformity=:L2)
    M = TestFESpace(Γ,
      reffeₗ;
      conformity=:L2,
      dirichlet_tags=collect(5:8))
    Y = MultiFieldFESpace([V, Q, M])

    # Define trial FEspaces
    U = TrialFESpace(V)
    P = TrialFESpace(Q)
    L = TrialFESpace(M, p)
    X = MultiFieldFESpace([U, P, L])

    # FE formulation params
    τ = 1.0 # HDG stab parameter

    degree = 2 * (order + 1)
    dΩ = Measure(Ω, degree)
    n = get_cell_normal_vector(∂K)
    nₒ = get_cell_owner_normal_vector(∂K)
    d∂K = Measure(∂K, degree)
    a((uh, ph, lh), (vh, qh, mh)) = ∫(vh ⋅ uh - (∇ ⋅ vh) * ph - ∇(qh) ⋅ uh)dΩ +
                                    ∫((vh ⋅ n) * lh)d∂K +
                                    #∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                                    ∫(qh * (uh ⋅ n))d∂K +
                                    ∫(τ * qh * ph * (n ⋅ nₒ))d∂K -
                                    ∫(τ * qh * lh * (n ⋅ nₒ))d∂K +
                                    #∫(mh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                                    ∫(mh * (uh ⋅ n))d∂K +
                                    ∫(τ * mh * ph * (n ⋅ nₒ))d∂K -
                                    ∫(τ * mh * lh * (n ⋅ nₒ))d∂K
    l((vh, qh, mh)) = ∫(vh ⋅ f + qh * (∇ ⋅ u)) * dΩ
    op = HybridAffineFEOperator((u, v) -> (a(u, v), l(v)), X, Y, [1, 2], [3])
    op
  end

  function solve_darcy_lhdg(n,order)
    partition = (0, 1, 0, 1)
    cells = (n, n)
    model = simplexify(CartesianDiscreteModel(partition, cells))
    D = num_cell_dims(model)
    Ω = Triangulation(ReferenceFE{D}, model)

    degree = 2
    dΩ = Measure(Ω, degree)
    ∂K = GridapHybrid.Skeleton(model)
    d∂K = Measure(∂K, degree)

    op = build_darcy_lhdg(model, 0)

    order_M0 = 1
    reffe_M0 = ReferenceFE(lagrangian, Float64, order_M0; space=:Q)
    M0 = TestFESpace(Ω, reffe_M0; conformity=:H1, dirichlet_tags="boundary")

    a(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ
    A0 = assemble_matrix(a, M0, M0)

    x1 = zeros(length(op.skeleton_op.op.vector))
    num_iter=nonnested_mg_2_level_v_cycle!(x1, op, A0, M0, dΩ, d∂K; rtol=1.0e-6, maxiter=1000)

    xexact=-op.skeleton_op.op.matrix\op.skeleton_op.op.vector
    num_iter,norm(xexact-x1)
  end

  iters=Int[]
  for n in (5,10,15,20,30,40,50)
    num_iter,_=solve_darcy_lhdg(n,0)
    append!(iters,num_iter)
  end
  println(iters)

end
