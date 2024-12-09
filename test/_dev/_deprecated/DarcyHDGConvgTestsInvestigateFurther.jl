################################################################################################
# u + ∇p = 0 in Ω;   ∇ ⋅ u = f in Ω; u ⋅ n = g on Γ, ∫ p(x) dx = 0 in Ω                        #
################################################################################################
# HDG Spaces                                                                                   #
################################################################################################
# Vₕ  = { vₕ ∈ L²(τₕ)ᵈ: v|ₖ in Pᵣ(K) ∀ K ∈ τₕ },
# Qₕ  = { qₕ ∈ L²(τₕ) : q|ₖ in Pᵣ(K) ∀ K ∈ τₕ },
# Mₕ  = { ̂mₕ ∈ L²(Εₕ) : ̂mₕ|ₑ in Pᵣ(e) ∀ e ∈ Εₕ },
################################################################################################
# Notations (unless explicitly stated)                                                                                    # 
################################################################################################
# (⋅, ⋅) - L²(τₕ) inner product
# ⟨⋅, ⋅⟩ - L²(∂τₕ) inner product; ∂τₕ ≐ { ∂K: K ∈ τₕ }
# n - normal vector 
# τ - HDG Stabilization parameter
#################################################################################################
# HDG weak formulation                                                                          #
#################################################################################################
# For all (vₕ, qₕ, ̂mₕ) ∈ Vₕ × Qₕ × Mₕ, find (pₕ, uₕ, ̂uₕ) ∈ Vₕ × Qₕ × Mₕ such that               #
# ( uₕ, vₕ ) - ( pₕ, ∇ ⋅ vₕ ) - ⟨ ̂pₕ, vₕ ⋅ n ⟩ = 0,                                             #
#            - ( uₕ, ∇ qₕ )   - ⟨ ̂uₕ ⋅ n, qₕ ⟩ = ( f, qₕ ),                                     #
#                               ⟨ ̂uₕ ⋅ n, ̂mₕ ⟩ = ⟨ g, ̂mₕ ⟩, ( inner-product over Γ )                     
#                               ⟨ ̂uₕ ⋅ n, ̂mₕ ⟩ = 0,      
#                                   ( pₕ, 1 )  = 0,                                                                 
# where,                                                                                        #
# ̂uₕ ⋅ n = uₕ ⋅ n - τ ( pₕ - ̂pₕ )  on ∂τₕ                                                       #
#################################################################################################                    


module DarcyHDGConvgTests

using Test
using Gridap
using FillArrays
using Gridap.Geometry
using GridapHybrid
using Plots

p(x) = sin(2*π*x[1])*sin(2*π*x[2])*x[1]*(x[1]-1)*x[2]*(x[2]-1)   # ex 1: g = 0 
# p(x) = sin(2*π*x[1])*sin(2*π*x[2]) # ex 2: u ⋅ n = g   (To Do) # ex 2: g non-zero
u(x) = -∇(p)(x)
n = 4
# function solve_darcy_hdg((n,n),order)
    # Geometry
    partition = (0,1,0,1)
    cells = (n,n)
    model = CartesianDiscreteModel(partition,cells)
    D = num_cell_dims(model)
    Ω = Triangulation(ReferenceFE{D},model)
    Γ = Triangulation(ReferenceFE{D-1},model)
    ∂K = GridapHybrid.Skeleton(model)

    # Reference FEs
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order;space=:P)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
    # reffeₚ = ReferenceFE(lagrangian,Float64,order;space=:P)
    reffeₗ = ReferenceFE(lagrangian,Float64,order;space=:P)

    # Define test FESpaces
    V = TestFESpace(Ω  , reffeᵤ; conformity=:L2)
    Q = TestFESpace(Ω  , reffeₚ; conformity=:L2, constraint=:zeromean)
    M = TestFESpace(Γ, reffeₗ; conformity=:L2)   
    Y = MultiFieldFESpace([V,Q,M])

    # Define trial FEspaces
    U = TrialFESpace(V)
    P = TrialFESpace(Q)
    L = TrialFESpace(M)
    X = MultiFieldFESpace([U, P, L])

    # FE formulation params
    τ = 1.0 # HDG stab parameter

    degree = 2*(order+1)
    dΩ     = Measure(Ω,degree)
    n      = get_cell_normal_vector(∂K)
    nₒ     = get_cell_owner_normal_vector(∂K)
    d∂K    = Measure(∂K,degree)

    yh = get_fe_basis(Y)
    xh = get_trial_fe_basis(X)

    (uh,ph,lh) = xh
    (vh,qh,mh) = yh

    f = (∇ ⋅ u)

    # Does not work: Gives error while computing op
    # ∂Ω = BoundaryTriangulation(model)
    # d∂Ω= Measure(∂Ω, degree)
    # nbd= get_normal_vector(∂Ω) 
    # g = u ⋅ nbd
    # op = ∫(g)*mh*d∂Ω

    # a((uh,ph,lh),(vh,qh,mh)) =∫( vh⋅uh - (∇⋅vh)*ph - ∇(qh)⋅uh )dΩ +    
    #                           ∫((vh⋅n)*lh)d∂K +                        
    #                           ∫(qh*(uh⋅n))d∂K -                        
    #                           ∫(τ*qh*ph*(n⋅nₒ))d∂K +                   
    #                           ∫(τ*qh*lh*(n⋅nₒ))d∂K +                   
    #                           ∫(mh*(uh⋅n))d∂K -
    #                           ∫(τ*mh*ph*(n⋅nₒ))d∂K +
    #                           ∫(τ*mh*lh*(n⋅nₒ))d∂K

    a((uh,ph,lh),(vh,qh,mh)) =∫( vh⋅uh - (∇⋅vh)*ph - ∇(qh)⋅uh )dΩ +    # This form gives a MatrixBlock{Matrix{VectorValue{2, Float64}}}
                              ∫((vh⋅n)*lh)d∂K +                        # compared to the above commented out form.
                              ∫(qh*(uh⋅n))d∂K -                        
                              ∫(τ*qh*ph)d∂K +                   
                              ∫(τ*qh*lh)d∂K +                   
                              ∫(mh*(uh⋅n))d∂K -
                              ∫(τ*mh*ph)d∂K +
                              ∫(τ*mh*lh)d∂K

    l((vh,qh,mh)) = ∫(qh*f)*dΩ

    op=HybridAffineFEOperator((u,v)->(a(u,v),l(v)), X, Y, [1,2], [3])
    xh=solve(op)

    uh, ph, _ = xh
    e = p - ph
    
  #   return sqrt(sum(∫(e⋅e)dΩ))

  # end
  
  function conv_test(ns,order)
    el2 = Float64[]
    hs = Float64[]
    for n in ns
      l2 = solve_darcy_hdg((n,n),order)
      println(l2)
      h = 1.0/n
      push!(el2,l2)
      push!(hs,h)
    end
    println(el2)
    return el2, hs
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x), 1), x) \ y
    linreg[2]
  end

  # ns=[8,16,32,64,128]
  ns=[8,16,32]
  order=1
  # el, hs = conv_test(ns,order)
  # println("Slope L2-norm u: $(slope(hs,el))")
  # slopek  =[Float64(ni)^(-(order)) for ni in ns]
  # display(plot(hs,[el slopek],
  #   xaxis=:log, yaxis=:log,
  #   label=["L2u (measured)" "slope k"],
  #   shape=:auto,
  #   xlabel="h",ylabel="L2 error",legend=:bottomright))

end # module