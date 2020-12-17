using Gridap

# As a first step, let us implement a so-called hybrid primal finite element method (see [Brezzi and Fortin]). This is equivalent to the Crouzieaux-Raviart FE space on triangular meshes (see [Ern and Guermond]).  Usually, the Crouzieaux-Raviart method is not implemented this way, but this is useful for us, because we will gain experience with hybrid formulations.E.g., we will define the local spaces using sub-meshes, in order to implement hybrid multiscale methods (e.g., analysed in [Barrenechea et al]), and be a little bit closer to a BDDC preconditioner.

# Manufactured solution
u(x) = x[1] + x[2]
f(x) = -Δ(u)

# Model
domain = (0,1,0,1)
cells = (3,3)
D = length(cells)
model_quad = CartesianDiscreteModel(domain,cells)
model = simplexify(model_quad)
model_facets = DiscreteModel(Polytope{D-1},model)

# FESpaces
k = 1
reffe_u = ReferenceFE(lagrangian,Float64,k)
reffe_λ = ReferenceFE(lagrangian,Float64,k-1)

V = FESpace(model,reffe_u,conformity=:L2)

# This function is not yet implemented. I think it should be possible to have L2 conformity and strong imposition of Dirichlet data, for `FESpace`s that are nodal-based, and thus allow one to define the concept of ownership of DOFs by geometrical entities
# M = TestFESpace(model_facets,reffe_λ,conformity=:L2,dirichlet_tags="boundary")

# In any case, for Dirichlet boundary conditions, this is not important. Since λ is a flux, i.e., n⋅∇u, the Dirichlet boundary conditions on the flux are the Neumann boundaries of the primal problem. So, for the moment, we leave this discussion aside.

Λ = TestFESpace(model_facets,reffe_λ,conformity=:L2)

# Integration machinery
Ω = Triangulation(model)
Γ = SkeletonTriangulation(model)
∂Ω = BoundaryTriangulation(model,tags="boundary")

nΓ = get_normal_vector(Γ)
n∂Ω = get_normal_vector(∂Ω)

degree = 2*k
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
d∂Ω = Measure(Γ,degree)

# In the following, we want to write the global problem. Later on, we will discuss how to eliminate (part of) the u unknown using local problems.

function a(x,y)
  u, λ = x
  v, μ = y
  # Assuming that the jump does not include the normal vector (?)
  ∫(∇(u)⊙∇(v))dT + ∫(jump(v) * λ)dΓ + ∫(jump(u) * μ)dΓ +
  + ∫(v * λ)d∂Ω + ∫(u * μ)d∂Ω
end

function b(y)
  v, μ = y
  ∫(v*f)*dT + ∫(μ*u)*d∂Ω
end

M = MultiFieldFESpace([V,Λ])
op = AffineFEOperator(a,l,M,M)

u, λ = solve(op)

E = Triangulation(model_facets)
writevtk(E,"lambda",cellfields=["λ"=>λ])
writevtk(T,"results",cellfields=["u"=>u])

# We can rewrite the same problem using a more cell-wise approach.
# Development 2: Define the `CellBoundaries` of a given mesh

∂T = CellBoundaries(model)
# or
# ∂T = CellBoundaries(Ω)
d∂T = Measure(∂T)
n∂T = get_normal_vector(∂T)

# With this, we could now implement the same method above but writing the form as follows.

function a_alt(u,λ)
  s = n∂T⋅n.⁺ # One way to define a +/- sign to each side of the facet
  ∫(∇(u)⊙∇(v))dT + ∫(s*v*λ)dΓ + ∫(s*u*μ)dΓ +
  + ∫(s*v*λ)d∂Ω + ∫(s*u*μ)d∂Ω
end

l_alt(v,μ) = ∫(v*f)dT + ∫(μ*u)d∂T

# Can we make it work? Do we get the same result?

# Just to check things, the previous space is the trace of the so-called Raviart-Thomas space (see [Brezzi and Fortin]). So, we can in fact implement the same problem without the need to define traces spaces on the skeleton of the mesh. This fact could be used to check that everything works as expected.

# reffe_q = ReferenceFE(raviart_thomas,Float64,k)

# Q = TestFESpace(model,reffe_q,conformity=:HDiv)

# # If we want to use the RT space, the bilinear forms would be the following.

# function a_rt(x,y)
#   u, q = x
#   v, r = y
#   # Assuming that the jump does not include the normal vector (?)
#   # I think that the normal being taken is the same (it just changes the sign of the flux)
#   ∫(∇(u)⊙∇(v))dT + ∫(jump(v) * q⋅nΓ.⁺)dΓ + ∫(jump(u) * q⋅nΓ.⁺)dΓ +
#   + ∫(v * q⋅n∂Ω)d∂Ω + ∫(u * q⋅n∂Ω)d∂Ω
# end

# function b_rt(y)
#   v, r = y
#   ∫(v*f)*dT + ∫(u * r⋅n∂Ω)*d∂Ω
# end

# SECOND STEP

# In any case, we would like to be able to write this problem eliminating locally most of the u DOFs, i.e., to use static condensation of local variables and implement other hybrid methods in the future. We will need something in the following lines.

# We can eliminate the interior DOFs at each cell, i.e., only
# put the constant of the primal variable and the flux in the global system.
# Check the details e.g. in [Barrenechea et al] or [Araya et al, eqs. (1.7) and (1.8)]. Even though these articles involve multiscale spaces, the formulation is independent of this. This implementation is slightly more involved, because it involves 1) a space V₀⟂ identical to V but without the constant term. That should not be that hard, e.g., instead of using the standard Lagrangian FEs, one can use the monomial basis without the 1 term. Assuming we have this...

# Now, we can try to implement the previous problem using the approach in [Araya et al]

reffe_u_0 = ReferenceFE(lagrangian,Float64,0)
V₀ = FESpace(model,reffe_u_0,conformity=:L2)

# We would need the space with zero mean value, not implemented
reffe_u_0_⟂ = ReferenceFE(zero_mean_lagrangian,Float64,k) # to be implemented
V₀⟂ = FESpace(model,reffe_u_0_⟂)

# Local problem 1
a_local_⟂(u,v) = ∫(∇(u)⊙∇(v))dT

function l_local_⟂_λ(λ)
  s = n∂T⋅n.⁺ # One way to define a +/- sign to each side of the facet
  l_local
  l(v) = -∫(s*v*λ)d∂Ω # Given λ, return a linear form for v
end

l_local_⟂_f(v) = ∫(v*f)dT

# Now we create a Map `LocalSolve` that computes (lazily) the solution of the following cell-wise problem (no global problem!)
# u ∈ V0⟂ : a_local_⟂(u,v) = l_local_⟂(v) ∀ v ∈ V0⟂
Tf = LocalSolve(a_local_⟂,l_local_⟂,V₀⟂,V₀⟂)
T∂(λ)) = LocalSolve(a_local_⟂,l_local_⟂_λ(λ),V₀⟂,V₀⟂)
# see the problems (1.7) and (1.8) in [Araya et al]

# The result should be a `CellField` in both cases, the second one after being evaluated at a given λ

# With these results, we can now implement

function a_2(x,y)
  u, λ = x
  v, μ = y
  # Assuming that the jump does not include the normal vector (?)
  ∫(jump(v) * λ)dΓ + ∫(jump(u+T∂(λ)) * μ)dΓ
end

function l_2(y)
  v, μ = y
  ∫(v*f)*dΩ + ∫(μ*u)*d∂Ω - ∫(μ*jump(Tf))*dΓ
  # Here, I see a problem, integrating this way we are computing Tf many times
end

M₀ = MultiFieldFESpace([V₀,Λ])
op = AffineFEOperator(a_2,l_2,M₀,M₀)

u₀, λ = solve(op)
u = u₀ + Tf + T∂(λ)

writevtk(E,"lambda",cellfields=["λ"=>λ])
writevtk(T,"results",cellfields=["u"=>u])

# REFERENCES

# [Barrenechea et al] https://doi.org/10.1007/s00211-020-01103-5
# [Araya et al] https://doi.org/10.1137/120888223
