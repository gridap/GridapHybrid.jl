
using Gridap
using GridapHybrid
using FillArrays

include("P_m.jl")

# Q0: M1->M0
# uh1 belongs to M1
function M1_to_M0(uh1, M0, dΩ, d∂K)
   M1=Gridap.FESpaces.get_fe_space(uh1)
   vh0=get_fe_basis(M0)
   uh0=get_trial_fe_basis(M0)
   dcM = ∫(vh0*uh0)d∂K
   M = assemble_matrix(dcM,M0,M0)

   ∂K_diam_array = get_array(∫(1)d∂K)
   K_diam_array  = get_array(∫(1)dΩ)
   ∂K_diam_cf    = CellField(∂K_diam_array, d∂K.quad.trian, ReferenceDomain())
   K_diam_cf     = CellField(K_diam_array , d∂K.quad.trian, ReferenceDomain())

   uh0_in_M1=M0_to_M1(vh0,M0,M1,d∂K)
   dcb = ∫(uh1*uh0_in_M1*(K_diam_cf/∂K_diam_cf))d∂K

   b = assemble_vector(dcb,M0)
   FEFunction(M0,M\b)
end

function _to_trial_basis(a::Gridap.FESpaces.FEBasis,M::FESpace)
  if Gridap.FESpaces.BasisStyle(a) == Gridap.FESpaces.TrialBasis()
    return a
  else
    get_trial_fe_basis(M)
  end
end

function _to_trial_basis(a::Gridap.FESpaces.FEFunction,M::FESpace)
  a
end

# I1: M0->M1
# uh0 can be either a FE function in M0 or a basis of the trial/test space
function M0_to_M1(uh0::Union{<:Gridap.FESpaces.FEBasis,<:Gridap.FESpaces.FEFunction},M0,M1,d∂K)
  T0=get_triangulation(M0)
  T1=get_triangulation(M1)
  model=get_background_model(T0)
  @assert model === get_background_model(T1)
  D=num_cell_dims(model)
  @assert isa(T0,Triangulation{D,D})
  @assert isa(T1,Triangulation{D-1,D})

  # If uh0 is a FE function/trial basis, leave it untouched.
  # If uh0 is a test basis, transform it to trial. Otherwise, the
  # RHS in the following system does not have the appropriate shape.
  uh0t = _to_trial_basis(uh0,M0)

  μ = get_fe_basis(M1)
  v = get_trial_fe_basis(M1)

  A = ∫(μ*v)d∂K
  B = ∫(μ*uh0t)d∂K

  # [c][f][f,f]
  A_array = _remove_sum_facets(_remove_densify(Gridap.CellData.get_array(A)))

  # Assuming uh0t is a FE basis (always will be trial)
  # [c][f][f]
  B_array = _remove_sum_facets(_remove_densify(Gridap.CellData.get_array(B)))

  # Assuming uh0t is a FE basis (always will be trial)
  # [c][f][1]
  pm_uh_dofs = lazy_map(GridapHybrid.compute_bulk_to_skeleton_l2_projection_dofs, A_array, B_array)

  # [c][f][f]
  if (isa(uh0,FEFunction))
    M1_basis=v
  else
    @assert isa(uh0,Gridap.FESpaces.FEBasis)
    if (Gridap.FESpaces.BasisStyle(uh0)==Gridap.FESpaces.TrialBasis())
      M1_basis=v
    else
      M1_basis=μ
    end
  end

  M1_basis_d∂K      = Gridap.CellData.change_domain(M1_basis, d∂K.quad.trian, ReferenceDomain())
  M1_basis_d∂K_data = Gridap.CellData.get_data(M1_basis_d∂K)

  # Assuming uh0t is a FE basis
  # [c][f][1]   pm_uh_dofs
  # [c][f][f]   v_d∂K_data (test)
  # [c][f][1,f] v_d∂K_data (trial)

  # [c][f][1]
  field_array =
    lazy_map(GridapHybrid.setup_bulk_to_skeleton_l2_projected_fields, pm_uh_dofs, M1_basis_d∂K_data)

  cf=Gridap.CellData.GenericCellField(field_array, d∂K.quad.trian, ReferenceDomain())
end


partition = (0,1,0,1)
cells     = (2,2)
model     = CartesianDiscreteModel(partition, cells)
D         = num_cell_dims(model)
Ω         = Triangulation(ReferenceFE{D}, model)
dΩ        = Measure(Ω,2)
Γ         = Triangulation(ReferenceFE{D-1}, model)
∂K        = GridapHybrid.Skeleton(model)
d∂K       = Measure(∂K,2)


order_M0 = 1
reffe_M0 = ReferenceFE(lagrangian, Float64, order_M0; space = :Q)

order_M1 = 0
reffe_M1 = ReferenceFE(lagrangian, Float64, order_M1; space = :Q)

M0 = TestFESpace(Ω, reffe_M0; conformity = :H1)
M1 = TestFESpace(Γ, reffe_M1; conformity = :L2)

uh0=FEFunction(M0,rand(num_free_dofs(M0)))
cf=M0_to_M1(uh0,M0,M1,d∂K)

uh0=get_fe_basis(M0)
cf=M0_to_M1(uh0,M0,M1,d∂K)

uh1=FEFunction(M1,ones(num_free_dofs(M1)))
uh0=M1_to_M0(uh1, M0, dΩ, d∂K)

writevtk(Ω, "kk", cellfields=["uh0"=>uh0])
