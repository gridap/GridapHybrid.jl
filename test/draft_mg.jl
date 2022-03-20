
using Gridap
using GridapHybrid
using FillArrays
using IterativeSolvers
using Printf

include("P_m.jl")

# Q0: M1->M0
# uh1 belongs to M1
function M1_to_M0(uh1, M0, dΩ, d∂K)
  M1 = Gridap.FESpaces.get_fe_space(uh1)
  vh0 = get_fe_basis(M0)
  uh0 = get_trial_fe_basis(M0)
  dcM = ∫(vh0 * uh0)d∂K
  M = assemble_matrix(dcM, M0, M0)

  ∂K_diam_array = get_array(∫(1)d∂K)
  K_diam_array = get_array(∫(1)dΩ)
  ∂K_diam_cf = CellField(∂K_diam_array, d∂K.quad.trian, ReferenceDomain())
  K_diam_cf = CellField(K_diam_array, d∂K.quad.trian, ReferenceDomain())

  uh0_in_M1 = M0_to_M1(vh0, M0, M1, d∂K)
  # dcb = ∫(uh1 * uh0_in_M1 * (K_diam_cf / ∂K_diam_cf))d∂K
  dcb = ∫(uh1 * uh0_in_M1)d∂K

  b = assemble_vector(dcb, M0)
  FEFunction(M0, M \ b)
end

function _to_trial_basis(a::Gridap.FESpaces.FEBasis, M::FESpace)
  if Gridap.FESpaces.BasisStyle(a) == Gridap.FESpaces.TrialBasis()
    return a
  else
    get_trial_fe_basis(M)
  end
end

function _to_trial_basis(a::Gridap.FESpaces.FEFunction, M::FESpace)
  a
end

# I1: M0->M1
# uh0 can be either a FE function in M0 or a basis of the trial/test space
function M0_to_M1(uh0::Union{<:Gridap.FESpaces.FEBasis,<:Gridap.FESpaces.FEFunction}, M0, M1, d∂K)
  T0 = get_triangulation(M0)
  T1 = get_triangulation(M1)
  model = get_background_model(T0)
  @assert model === get_background_model(T1)
  D = num_cell_dims(model)
  @assert isa(T0, Triangulation{D,D})
  @assert isa(T1, Triangulation{D - 1,D})

  # If uh0 is a FE function/trial basis, leave it untouched.
  # If uh0 is a test basis, transform it to trial. Otherwise, the
  # RHS in the following system does not have the appropriate shape.
  uh0t = _to_trial_basis(uh0, M0)

  μ = get_fe_basis(M1)
  v = get_trial_fe_basis(M1)

  A = ∫(μ * v)d∂K
  B = ∫(μ * uh0t)d∂K

  # [c][f][f,f]
  A_array = _remove_sum_facets(_remove_densify(Gridap.CellData.get_array(A)))

  # Assuming uh0t is a FE basis (always will be trial)
  # [c][f][f]
  B_array = _remove_sum_facets(_remove_densify(Gridap.CellData.get_array(B)))

  # Assuming uh0t is a FE basis (always will be trial)
  # [c][f][1]
  pm_uh_dofs = lazy_map(GridapHybrid.compute_bulk_to_skeleton_l2_projection_dofs, A_array, B_array)

  # [c][f][f]
  if (isa(uh0, FEFunction))
    M1_basis = v
  else
    @assert isa(uh0, Gridap.FESpaces.FEBasis)
    if (Gridap.FESpaces.BasisStyle(uh0) == Gridap.FESpaces.TrialBasis())
      M1_basis = v
    else
      M1_basis = μ
    end
  end

  M1_basis_d∂K = Gridap.CellData.change_domain(M1_basis, d∂K.quad.trian, ReferenceDomain())
  M1_basis_d∂K_data = Gridap.CellData.get_data(M1_basis_d∂K)

  # Assuming uh0t is a FE basis
  # [c][f][1]   pm_uh_dofs
  # [c][f][f]   v_d∂K_data (test)
  # [c][f][1,f] v_d∂K_data (trial)

  # [c][f][1]
  field_array =
    lazy_map(GridapHybrid.setup_bulk_to_skeleton_l2_projected_fields, pm_uh_dofs, M1_basis_d∂K_data)

  cf = Gridap.CellData.GenericCellField(field_array, d∂K.quad.trian, ReferenceDomain())
end

# op  : Matrix + RHS  corresponding to the HDG skeleton problem
#       (after static condensation of bulk unknowns)
# A0  : operator in coarse space
# M0  : coarse space
# dΩ  : bulk measure
# d∂K : skeleton measure
function nonnested_mg_2_level_v_cycle!(x, op1, A0, M0, dΩ, d∂K; rtol=1.0e-06, maxiter=10)
  # TO-DO: - sign should not be here.
  #        for unknown reason, matrix is symmetric negative definite
  A1 = -op1.skeleton_op.op.matrix
  b1 = op1.skeleton_op.op.vector

  M1 = op1.skeleton_op.test

  r = b1 - A1 * x  # dynamic memory alloc
  nrm_r0 = norm(r)
  nrm_r = nrm_r0
  current_iter = 0
  rel_res = nrm_r / nrm_r0
  @printf "%6s  %12s" "Iter" "Rel res\n"
  while current_iter <= maxiter && rel_res > rtol
    @printf "%6i  %12.4e\n" current_iter rel_res
    smooth!(x, A1, b1; maxiter=2)
    r .= b1 .- A1 * x  # dynamic memory alloc
    I1_B0_Q0!(x, M1, A1, r, A0, M0, dΩ, d∂K)
    smooth!(x, A1, b1; maxiter=2)
    r .= b1 .- A1 * x
    nrm_r = norm(r)
    rel_res = nrm_r / nrm_r0
    current_iter += 1
  end
  return current_iter
end

function I1_B0_Q0!(x, M1, A1, b1, A0, M0, dΩ, d∂K)
  uh1 = FEFunction(M1, b1)
  uh0 = M1_to_M0(uh1, M0, dΩ, d∂K)
  x0 = A0 \ Gridap.FESpaces.get_free_dof_values(uh0)
  uh0 = FEFunction(M0, x0)
  δh1_∂K = M0_to_M1(uh0, M0, M1, d∂K)

  # To-Do: In the following bunch of lines I am doing something
  # numerically that can be just be performed symbollicaly.
  # Namely to transform a SkeletonArray to a FacetArray.
  a(u, v) = ∫(v * u)d∂K
  l(v) = ∫(v * δh1_∂K)d∂K
  Aδh1 = assemble_matrix(a, M1, M1)
  bδh1 = assemble_vector(l, M1)
  δh1 = FEFunction(M1, Aδh1 \ bδh1)

  x .= x .+ Gridap.FESpaces.get_free_dof_values(δh1)
end

function smooth!(x, A, b; maxiter=1)
  gauss_seidel!(x, A, b; maxiter=maxiter)
end


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
