using Gridap
using GridapHybrid

function Gridap.Arrays.return_cache(::typeof(\),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  c = a\b
  Gridap.Arrays.CachedArray(c)
end

function Gridap.Arrays.evaluate!(cache,::typeof(\),a::AbstractMatrix{<:Number},b::AbstractVector{<:Number})
  m = size(a,1)
  Gridap.Arrays.setsize!(cache,(m,))
  c = cache.array
  c .= a\b
  c
end

function Gridap.Arrays.evaluate!(cache,::typeof(\),a::AbstractMatrix{<:Number},b::AbstractMatrix{<:Number})
  m = size(a,1)
  n = size(b,2)
  Gridap.Arrays.setsize!(cache,(m,n))
  c = cache.array
  c .= a\b
  c
end


# TO-THINK: Perhaps a better name??? E.g., Local FE problem???
struct LocalAffineFEOperator{T1,T2,T3,T4} <: GridapType
  LHS_form               :: T1
  RHS_form               :: T2
  trial_space            :: T3
  test_space             :: T4
  LHS_contribs           :: Gridap.CellData.DomainContribution

  function LocalAffineFEOperator(weakform,
                               trial_space :: A,
                               test_space  :: B) where {A,B}

     us = get_trial_fe_basis(trial_space)
     vs = get_fe_basis(test_space)

     LHS_form, RHS_form = weakform
     LHS_contribs       = LHS_form(us,vs)

     T1=typeof(LHS_form)
     T2=typeof(RHS_form)
     T3=typeof(trial_space)
     T4=typeof(test_space)
     new{T1,T2,T3,T4}(LHS_form,
                      RHS_form,
                      trial_space,
                      test_space,
                      LHS_contribs)
  end
end

# TO-THINK: MultiFieldFEBasis are NOT of type <: FEBasis.
#           They are of type MultiFieldCellField. Besides,
#           I have the feeling the code below only works for
#           SingleFieldFEBasis.

function (::LocalAffineFEOperator)(bases::Gridap.FESpaces.FEBasis...)
  basis_style=Gridap.FESpaces.BasisStyle(first(bases))
  Gridap.Helpers.@check all( [basis_style==Gridap.FESpaces.BasisStyle(basis) for basis in bases ] )

  lhs_skeleton=GridapHybrid._find_skeleton(op.LHS_contribs)
  if length(lhs_skeleton)==1
    LHS_contribs = Gridap.Hybrid._merge_bulk_and_skeleton_contributions(op.LHS_contribs)
  else
    Gridap.Helpers.@check length(lhs_skeleton)==0
    LHS_contribs = op.LHS_contribs
  end

  vs = get_fe_basis(op.test_space)
  RHS_contribs = op.RHS_form(bases,vs)

  rhs_skeleton=GridapHybrid._find_skeleton(RHS_contribs)
  if length(rhs_skeleton)==1
    RHS_contribs = GridapHybrid._merge_bulk_and_skeleton_matrix_contributions(RHS_contribs)
  end
  A_array = get_array(LHS_contribs)
  B_array = get_array(RHS_contribs)
  dofs = lazy_map(\,A_array,B_array)

  # We get the test space so that we dont deal with
  # the transpose underlying the trial basis
  basis=get_fe_basis(op.test_space)

  range_basis_data = lazy_map(Gridap.Fields.linear_combination,
                              dofs,
                              Gridap.CellData.get_data(basis))

  if (basis_style==Gridap.FESpaces.TrialBasis())
    range_basis_data = lazy_map(transpose,range_basis_data)
  end

  Gridap.FESpaces.similar_fe_basis(basis,
                                   range_basis_data,
                                   Gridap.CellData.get_domains(RHS_contribs)...,
                                   basis_style,
                                   DomainStyle(basis))
end


model=CartesianDiscreteModel((0,1,0,1),(2,2))

# Geometry
D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D},model)
Γ = Triangulation(ReferenceFE{D-1},model)
∂K = GridapHybrid.Skeleton(model)

# Reference FEs
order  = 1
refferecᵤ = ReferenceFE(lagrangian,Float64,order+1;space=:P)
reffeᵤ    = ReferenceFE(lagrangian,Float64,order  ;space=:P)

# Define test FESpaces
VKR  = TestFESpace(Ω  , refferecᵤ; conformity=:L2)
VK   = TestFESpace(Ω  , reffeᵤ   ; conformity=:L2)
V∂K  = TestFESpace(Γ  , reffeᵤ   ; conformity=:L2)


# Define trial FEspaces
UKR = TrialFESpace(VKR)
UK  = TrialFESpace(VK)
U∂K = TrialFESpace(V∂K)

degree = 2*order+1
dΩ     = Measure(Ω,degree)
nK     = get_cell_normal_vector(∂K)
d∂K    = Measure(∂K,degree)

v=get_fe_basis(VKR)
u∂K=get_trial_fe_basis(U∂K)
dc=∫((∇(v)⋅nK)*u∂K)d∂K


a(u,v)=∫(∇(v)⋅∇(u))dΩ
l( (uK,u∂K), v)=∫((∇(v)⋅nK)*u∂K)d∂K+∫(-Δ(v)*uK)dΩ

op=LocalAffineFEOperator((a,l),UKR,VKR)
UKR_basis=op(get_trial_fe_basis.((UK,U∂K))...)
