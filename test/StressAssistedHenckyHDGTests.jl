module StressAssistedHenckyHDGTests

# Implements HDG formulation for stress-assisted diffusion
# coupled problem with Hencky elasticity

using Test
using Gridap
using FillArrays
using GridapHybrid
using Plots

#const ν = 0.3 # Poisson ratio
#const E = 1.0 # Young's modulus
# Lame parameters
const λ = 1.0e+04 #(E*ν)/((1.0+ν)*(1.0-2.0*ν))  #
const μ = 1.0     #E/(2.0*(1.0+ν)) #

# Nonlinear constitutive form.
# ρ is the modulus-squared of the deviatoric strain tensor
function ζ(ρ)
  0.75*(1.0/(sqrt(1.0+ρ)))
end

function u_exact(x)
  #u1=0.1*(-cos(x[1])*sin(x[2])+((0.5*x[1]*x[1])/λ))
  #u2=0.1*(sin(x[1])*cos(x[2])+((0.5*x[2]*x[2])/λ))
  u1=x[1]+x[2]
  u2=2*x[1]-x[1]
  VectorValue((u1,u2))
end

function t_exact(x)
  t=0.5*(∇(u_exact)(x) + transpose(∇(u_exact)(x)))
  D=size(t)[1]
  vals = Gridap.TensorValues._flatten_upper_triangle(reshape(collect(t.data),(D,D)), Val(D))
  Gridap.TensorValues.SymTensorValue(vals)
end

function ϕ_exact(x)
  #cos(pi*x[1])*cos(pi*x[2])
  x[1]
end

function ω_exact(x)
  ∇(ϕ_exact)(x)
end

# To define nonlinear constitutive function ...
# σ: stress (tensor field)
# ϕ: solute concentration (scalar field)
function κ(σ,ϕ)
  #1.0+exp(tr(σ))+exp(-ϕ)
  one(σ)
end

# To define nonlinear constitutive function ...
# t: strain (tensor field)
function N(t)
  # td   = _dev(t)
  # tdsq = td⊙td
  # I    = one(t)
  # (λ + 2.0/3.0*ζ(tdsq))*tr(t)*I + 2.0*μ*(1.0-ζ(tdsq))*t
  t
end

# Extracts deviatoric component of t
function _dev(t)
  I = one(t)
  D = size(t)[1]
  invD=1.0/D
  t-(invD*tr(t))*I
end


C0=1.0
C1=1.0
C2=2.0
C3=1.0
# Solute dependent active stress component
function g(ϕ)
  # C0+C0*(ϕ^C1/(C2+ϕ^C3))
  # C0+(ϕ/(C2+ϕ))
  ϕ
end

# Solute concentration Dirichlet boundary
function ϕ₀(x)
  ϕ_exact(x)
end

# Displacement Dirichlet boundary
function u₀(x)
  u_exact(x)
end

function σ_exact(x)
  t=0.5*(∇(u_exact)(x) + transpose(∇(u_exact)(x)))
  res=N(t)-g(ϕ_exact(x))*one(t)
  D=size(res)[1]
  vals = Gridap.TensorValues._flatten_upper_triangle(reshape(collect(res.data),(D,D)), Val(D))
  Gridap.TensorValues.SymTensorValue(vals)
end

function ρ_exact(x)
  κ(σ_exact(x),ϕ_exact(x))⋅ω_exact(x)
end


# Vector field of body loads
function f(x)
  -(∇⋅(σ_exact))(x)
end

# Source/sink of solute concentration
function ℓ(x)
  0.0
end

# Appears in RHS of concentration gradient equation
function α(x)
  1.0
end

include("P_m.jl")

function solve_stress_assisted_diffusion_hencky_hdg(cells,order;write_results=false)
    # Geometry
    partition = (0,1,0,1)
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
    degree = 4 * (order + 1) # TO-DO: To think which is the minimum degree required
    dΩ = Measure(Ω, degree)
    n = get_cell_normal_vector(∂K)
    d∂K = Measure(∂K, degree)

    # Stabilization operator (solute concentration)
    τϕΓ = Float64(cells[1]) # cells[1] # To be defined???
    # Stabilization operator (displacements)
    τuΓ = Float64(cells[1]) # To be defined???

    println("h=$(1.0/cells[1]) τuΓ=$(τuΓ) τϕΓ=$(τϕΓ)")

    Y_basis    = get_fe_basis(Y)
    Y_TR_basis = get_trial_fe_basis(Y_TR)
    (_,_,_,_,_,_,_,uhΓ_basis) = Y_TR_basis

    # Testing for correctness of Pₘ projection
    xh              = FEFunction(Y_TR,rand(num_free_dofs(Y_TR)))
    _,_,_,_,_,_,_,μ = Y_basis
    _, uh, _, _ = xh
    dc=∫(μ⋅(uh-Pₘ(uh, uhΓ_basis, μ, d∂K)))d∂K
    for k in dc.dict[d∂K.quad.trian]
      @test all(k[8] .< 1.0e-14)
    end

    # Residual
    # current_iterate  => (ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ)
    # test space basis => (rh,nh,ψh,sh,τh,vh,ηh,μh)
    function residual((ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ),(rh,nh,ψh,sh,τh,vh,ηh,μh))
      Pm_uh=Pₘ(uh, uhΓ_basis, μh, d∂K)
      ∫( rh⋅((κ∘(σh,ϕh))⋅ωh) - rh⋅ρh )dΩ +
      ∫( nh⋅ωh+(∇⋅nh)*ϕh )dΩ - ∫((nh⋅n)*ϕhΓ)d∂K +
      ∫( ψh*(∇⋅ρh) - α*ψh*tr(th) - ψh*ℓ )dΩ +
           ∫( τϕΓ*ψh*ϕh )d∂K - ∫( τϕΓ*ψh*ϕhΓ )d∂K +
      ∫( sh⊙(N∘th) - sh⊙σh - tr(sh)*(g∘ϕh) )dΩ +
      ∫( τh⊙th + (∇⋅τh)⋅uh )dΩ - ∫((τh⋅n)⋅uhΓ)d∂K +
      ∫( ∇(vh)⊙σh - vh⋅f)dΩ - ∫(vh⋅(σh⋅n))d∂K + ∫(τuΓ*(vh⋅Pm_uh))d∂K -
                                                ∫(τuΓ*(vh⋅uhΓ))d∂K +
      ∫( ηh*(ωh⋅n) )d∂K + ∫( τϕΓ*ηh*ϕh )d∂K - ∫( τϕΓ*ηh*ϕhΓ )d∂K +
      ∫( μh⋅(σh⋅n) )d∂K - ∫( τuΓ*μh⋅Pm_uh)d∂K + ∫( τuΓ*μh⋅uhΓ )d∂K
    end

    # free_values = zeros(num_free_dofs(Y_TR))
    # xh          = FEFunction(Y_TR,free_values)
    # (ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ) = xh
    # (rh,nh,ψh,sh,τh,vh,ηh,μh)   = Y_basis
    # dc=residual((ωh,ρh,ϕh,th,σh,uh,ϕhΓ,uhΓ),(rh,nh,ψh,sh,τh,vh,ηh,μh))
    # res=assemble_vector(dc,Y)
    # dc=jacobian(x->residual(x,Y_basis),xh)
    # A=assemble_matrix(dc,Y_TR,Y)
    # op=HybridFEOperator(residual,Y_TR,Y,collect(1:6),[7,8])

    op=FEOperator(residual,Y_TR,Y)

    nls = NLSolver(show_trace=true, method=:newton)
    solver = FESolver(nls)

    xh0 = FEFunction(Y_TR,zeros(num_free_dofs(Y_TR)))
    xh, = solve!(xh0,solver,op)

    ωh,ρh,ϕh,th,σh,uh,_ = xh

    eω = ω_exact - ωh
    eρ = ρ_exact - ρh
    eϕ = ϕ_exact - ϕh
    et = t_exact - th
    eσ = σ_exact - σh
    eu = u_exact - uh
    norm2_ω=sqrt(sum(∫(eω ⊙ eω)dΩ))
    norm2_ρ=sqrt(sum(∫(eρ ⊙ eρ)dΩ))
    norm2_t=sqrt(sum(∫(et ⊙ et)dΩ))
    norm2_ϕ=sqrt(sum(∫(eϕ ⊙ eϕ)dΩ))
    norm2_σ=sqrt(sum(∫(eσ ⊙ eσ)dΩ))
    norm2_u=sqrt(sum(∫(eu ⋅ eu)dΩ))

    if (write_results)
       writevtk(Ω,"results_$(cells)_k=$(order)",
            cellfields=["σh"=>σh,"uh"=>uh,"ϕh"=>ϕh,
                        "σ_exact"=>σ_exact, "u_exact"=>u_exact, "ϕ_exact"=>ϕ_exact,
                        "eσ"=>eσ, "eu"=>eu,"eϕ"=>eϕ])
     end

     norm2_ω,norm2_ρ,norm2_t,norm2_ϕ,norm2_σ,norm2_u
    #@test sqrt(sum(∫(eu ⋅ eu)dΩ)) < 1.0e-12
    #@test sqrt(sum(∫(eσ ⊙ eσ)dΩ)) < 1.0e-12
  end

  function conv_test(ns,order)
    el2σ = Float64[]
    el2u = Float64[]
    hs = Float64[]
    for n in ns
      l2σ, l2u, _ = solve_stress_assisted_diffusion_hencky_hdg((n,n),order)
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

  ns=[8,16,32,64]
  order=1
  el2σ_PM, el2u_PM, hs = conv_test([8,16,32,64],1)
  plot(hs,[el2σ_PM el2u_PM slopek slopekp1 slopekp2],
    xaxis=:log, yaxis=:log,
    label=["L2σ (measured)" "L2u (measured)" "slope k" "slope k+1" "slope k+2"],
    shape=:auto,
    xlabel="h",ylabel="L2 error norm (PM)",legend=:bottomright)

  println("Slope L2-norm stress (PM): $(slope(hs,el2σ_PM))")
  println("Slope L2-norm      u (PM): $(slope(hs,el2u_PM))")

  solve_stress_assisted_diffusion_hencky_hdg((10,10),1,write_results=true)

end # module
