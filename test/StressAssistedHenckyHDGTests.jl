module StressAssistedHenckyHDGTests

# Implements HDG formulation for stress-assisted diffusion
# coupled problem with Hencky elasticity

using Test
using Gridap
using FillArrays
using Gridap.Geometry
using GridapHybrid
using Plots

include("tmp_fixes_symmetric_valued_lagrangian_reffes.jl")

const ν = 0.3 # Poisson ratio
const E = 1.0 # Young's modulus
# Lame parameters
const λ = (E*ν)/((1.0+ν)*(1.0-2.0*ν))  #
const μ = E/(2.0*(1.0+ν)) #


# Nonlinear constitutive form.
# ρ is the modulus-squared of the deviatoric stress tensor
function ζ(ρ)
  # ??? To be defined
  (3.0/4.0)*1.0/(sqrt(1.0+ρ))
end

# To define nonlinear constitutive function ...
# σ: stress (tensor field)
# ϕ: solute concentration (scalar field)
function κ(σ,ϕ)
  # ??? To be defined
  σ
end

# To define nonlinear constitutive function ...
# t: strain (tensor field)
function N(t)
  td   = _dev(t)
  tdsq = td⊙td
  I    = one(t)
  (λ + 2.0/3.0*ζ(tdsq))*tr(t)*I + 2.0*μ*(1.0-ζ(tdsq))*t
end

# Extracts deviatoric component of t
function _dev(t)
  I = one(t)
  D = size(t)[1]
  t-(tr(t)/D)*I
end

# Solute dependent active stress component
function g(ϕ)
  ϕ
  # ??? To be defined
end

# Solute concentration Dirichlet boundary
function ϕ₀(x)
  x[1]
  # ??? To be defined
end

# Displacement Dirichlet boundary
function u₀(x)
  x
  # ??? To be defined
end

# Vector field of body loads
function f(x)
  # (∇ ⋅ (σ_exact))(x)
  x
  # ??? To be defined
end

# Source/sink of solute concentration
function ℓ(x)
  x[1]
  # ??? To be defined
end

# ??? Appears in RHS of concentration gradient equation
function α(x)
  1.0
  # ??? To be defined
end

include("P_m.jl")

function Gridap.CellData.get_triangulation(f::Gridap.MultiField.MultiFieldCellField)
  s1 = first(f.single_fields)
  trian = get_triangulation(s1)
  trian
end

#function solve_stress_assisted_diffusion_hencky_hdg(cells,order;write_results=false)
    # Geometry
    partition = (0,1,0,1)
    cells=(2,2)
    order=1
    model = CartesianDiscreteModel(partition, cells)
    D = num_cell_dims(model)
    Ω = Triangulation(ReferenceFE{D}, model)
    Γ = Triangulation(ReferenceFE{D-1}, model)
    ∂K = GridapHybrid.Skeleton(model)

    # Reference FEs
    num_comp = D * (D + 1) ÷ 2 # Number of components of a symmetric tensor in D-dim space
    sym_tensor_type = Gridap.Fields.SymTensorValue{D,Float64,num_comp}

    # Concentration gradient space
    reffeω = ReferenceFE(lagrangian, VectorValue{D,Float64}, order; space = :P)
    # Concentration flux space
    reffeρ = reffeω
    # Solute concentration space
    reffeϕ = ReferenceFE(lagrangian, Float64, order; space = :P)
    # Trace of solute concentration space
    reffeϕΓ = reffeϕ
    # Trace of solute concentration space
    reffet  = ReferenceFE(lagrangian, sym_tensor_type, order; space = :P)
    # Symmetric stresses space
    reffeσ  = reffet
    # Displacement Space
    reffeu  = ReferenceFE(lagrangian, VectorValue{D,Float64}, order+1; space = :P)
    # Trace of Displacements Space
    reffeuΓ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order; space = :P)

    # Define test FESpaces
    W       = TestFESpace(Ω, reffeω; conformity = :L2)
    Ψ       = TestFESpace(Ω, reffeϕ; conformity = :L2)
    S       = TestFESpace(Ω, reffet; conformity = :L2)
    V       = TestFESpace(Ω, reffeu; conformity = :L2)
    MϕΓ     = TestFESpace(Γ, reffeϕΓ; conformity=:L2, dirichlet_tags=collect(5:8))
    MuΓ     = TestFESpace(Γ, reffeuΓ; conformity=:L2, dirichlet_tags=collect(5:8))
    Y       = MultiFieldFESpace([W, W, Ψ, S, S, V, MϕΓ, MuΓ])

    # Define trial FEspaces
    W_TR   = TrialFESpace(W)
    Ψ_TR   = TrialFESpace(Ψ)
    S_TR   = TrialFESpace(S)
    V_TR   = TrialFESpace(V)
    MϕΓ_TR = TrialFESpace(MϕΓ,ϕ₀)
    MuΓ_TR = TrialFESpace(MuΓ,u₀)
    Y_TR   = MultiFieldFESpace([W_TR, W_TR, Ψ_TR, S_TR, S_TR, V_TR, MϕΓ_TR, MuΓ_TR])

    # FE formulation params
    τ = cells[1] # alpha*order*order*(1.0/cells[1]) # HDG stab parameter
    println("h=$(1.0/cells[1]) τ=$(τ)")

    degree = 2 * (order + 1) # TO-DO: To think which is the minimum degree required
    dΩ = Measure(Ω, degree)
    n = get_cell_normal_vector(∂K)
    d∂K = Measure(∂K, degree)

    # Stabilization operator (solute concentration)
    τϕΓ = 1.0 # To be defined???
    # Stabilization operator (displacements)
    τuΓ = 1.0 # To be defined???

    Y_basis    = get_fe_basis(Y)
    Y_TR_basis = get_trial_fe_basis(Y_TR)
    (_,_,_,_,_,_,_,uhΓ_basis) = Y_TR_basis

    # Residual
    # current_iterate  => (ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ)
    # test space basis => (rh,nh,ψh,sh,τh,vh,ηh,μh)
    function residual((ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ),(rh,nh,ψh,sh,τh,vh,ηh,μh))
      ∫( rh⋅((κ∘(σh,ϕh))⋅ωh) - rh⋅ρh )dΩ +
      ∫( nh⋅ωh-(∇⋅nh)*ϕh )dΩ + ∫((nh⋅n)*ϕhΓ)d∂K  +
      ∫( ψh*(∇⋅ρh) - α*ψh*tr(th) - ψh*ℓ )dΩ +
      ∫( sh⊙(N∘th) - sh⊙σh - tr(sh)*(g∘ϕh) )dΩ +
      ∫( τh⊙th - (∇⋅τh)⋅uh )dΩ - ∫((τh⋅n)⋅uhΓ)d∂K +
      ∫( ∇(vh)⊙σh - vh⋅f)dΩ - ∫(vh⋅(σh⋅n))d∂K + ∫(τuΓ*(vh⋅Pₘ(uh, uhΓ_basis, μh, d∂K)))d∂K -
                                                ∫(τuΓ*(vh⋅uhΓ))d∂K +
      ∫( ηh*(ωh⋅n) )d∂K + ∫( τϕΓ*ηh*ϕh )d∂K - ∫( τϕΓ*ηh*ϕhΓ )d∂K +
      ∫( μh⋅(σh⋅n) )d∂K - ∫( τuΓ*μh⋅Pₘ(uh, uhΓ_basis, μh, d∂K) )d∂K + ∫( τuΓ*μh⋅uhΓ )d∂K
    end

    free_values = rand(num_free_dofs(Y_TR))
    xh          = FEFunction(Y_TR,free_values)

    (ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ) = xh
    (rh,nh,ψh,sh,τh,vh,ηh,μh)   = Y_basis
    dc=residual((ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ),(rh,nh,ψh,sh,τh,vh,ηh,μh))

    res=assemble_vector(dc,Y)

    # Jacobian
    # We aim at computing it using Automatic Differentiation
    # TO-DO: Need to think how to implement HybridFEOperator
    # (the counterpart of FEOperator for hybridizable nonlinear systems)

    # Bulk     unknowns -> [1,2,3,4,5,6]
    # Skeleton unknowns -> [7,8]

    # op = HybridAffineFEOperator((u, v) -> (a(u, v), l(v)), X, Y, [1, 2], [3])
    # xh = solve(op)

    # σh, uh, uhΓ = xh
    # eσ = σ_exact - σh
    # eu = u_exact - uh
    # norm2_σ=sqrt(sum(∫(eσ ⊙ eσ)dΩ))
    # norm2_u=sqrt(sum(∫(eu ⋅ eu)dΩ))
    # norm2_σ,norm2_u

    # if (write_results)
    #   writevtk(Ω,"results_$(cells)_k=$(order)",cellfields=["σh"=>σh,"uh"=>uh,"σ_exact"=>σ_exact, "error"=>eσ])
    # end
    # return     norm2_σ,norm2_u
    #@test sqrt(sum(∫(eu ⋅ eu)dΩ)) < 1.0e-12
    #@test sqrt(sum(∫(eσ ⊙ eσ)dΩ)) < 1.0e-12
  # end


  function conv_test(ns,order)
    el2σ = Float64[]
    el2u = Float64[]
    hs = Float64[]
    for n in ns
      l2σ, l2u = solve_stress_assisted_diffusion_hencky_hdg((n,n),order)
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

  cells=(2,2)
  order=1
  solve_stress_assisted_diffusion_hencky_hdg(cells,order;write_results=false)



  # ns=[8,16,32,64,128]
  # order=1
  # el2σ_PM, el2u_PM, hs = conv_test([8,16,32,64,128],1)
  # plot(hs,[el2σ_PM el2u_PM slopek slopekp1 slopekp2],
  #   xaxis=:log, yaxis=:log,
  #   label=["L2σ (measured)" "L2u (measured)" "slope k" "slope k+1" "slope k+2"],
  #   shape=:auto,
  #   xlabel="h",ylabel="L2 error norm (PM)",legend=:bottomright)

  # println("Slope L2-norm stress (PM): $(slope(hs,el2σ_PM))")
  # println("Slope L2-norm      u (PM): $(slope(hs,el2u_PM))")

end # module
