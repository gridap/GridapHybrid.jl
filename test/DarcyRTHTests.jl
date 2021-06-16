module IntegrationTests

using Test
using Gridap
using FillArrays
using Gridap.Geometry
using ExploringGridapHybridization

u(x) = VectorValue(1+x[1],1+x[2])
Gridap.divergence(::typeof(u)) = (x) -> 2
p(x) = -3.14
∇p(x) = VectorValue(0,0)
Gridap.∇(::typeof(p)) = ∇p
f(x) = u(x) + ∇p(x)
# Normal component of u(x) on Neumann boundary
function g(x)
  tol=1.0e-14
  if (abs(x[2])<tol)
    return -x[2] #-x[1]-x[2]
  elseif (abs(x[2]-1.0)<tol)
    return x[2] # x[1]+x[2]
  end
  Gridap.Helpers.@check false
end

function solve_darcy_rt_hdiv(model,order)
  V = FESpace(model,
              ReferenceFE(raviart_thomas,Float64,order),
              conformity=:Hdiv)#,dirichlet_tags=[5,6,7,8])
  Q = FESpace(model,ReferenceFE(lagrangian,Float64,order); conformity=:L2)
  U = TrialFESpace(V,u)
  P = TrialFESpace(Q)
  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])
  trian = Triangulation(model)
  degree = 2*(order+1)
  dΩ = Measure(trian,degree)
  neumanntags = [5,6,7,8]
  btrian = BoundaryTriangulation(model,tags=neumanntags)
  dΓ = Measure(btrian,degree)
  nb = get_normal_vector(btrian)
  a((u, p),(v, q)) = ∫( u⋅v - (∇⋅v)*p + q*(∇⋅u) )*dΩ
  b(( v, q)) = ∫( v⋅f + q*(∇⋅u))*dΩ - ∫((v⋅nb)*p )*dΓ
  op = AffineFEOperator(a,b,X,Y)
  #println(op.op.matrix)
  #println(op.op.vector)
  xh = solve(op)
end

function solve_darcy_hybrid_rt(model,order)

  # Geometry part
  D=2
  model_Γ = BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))

  # Functional part
  # To investigate what is needed to have an inf-sup stable triplet
  # for the RT-H method
  reffeᵤ = ReferenceFE(raviart_thomas,Float64,order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order)
  reffeₗ = ReferenceFE(lagrangian,Float64,order)

  # Compute the Dof values of Dirichlet DoFs (L2 projection)
  M = TestFESpace(model_Γ, reffeₗ; conformity=:L2, dirichlet_tags=[9])
  L = TrialFESpace(M)
  dirichlettags=[5,6,7,8]
  dirichlettrian=BoundaryTriangulation(model,tags=dirichlettags)
  degree = 2*(order+1)
  dΓd = Measure(dirichlettrian,degree)
  mh = get_fe_basis(M)
  lh = get_fe_basis(L)
  mass=∫(mh*lh)*dΓd
  rhs=∫(mh*p)*dΓd
  fdofsd=get_cell_dof_ids(M,dirichlettrian)
  data_a=Gridap.CellData.get_contribution(mass,dirichlettrian)
  data_b=Gridap.CellData.get_contribution(rhs,dirichlettrian)
  assem=SparseMatrixAssembler(M,L)
  A,b=assemble_matrix_and_vector(assem,(([],[],[]),
                                      ([data_a],[fdofsd],[fdofsd]),
                                      ([data_b],[fdofsd])))
  dirichlet_dofs=A\b

  M = TestFESpace(model_Γ, reffeₗ; conformity=:L2,dirichlet_tags=[5,6,7,8])
  fdofsd_new=get_cell_dof_ids(M,dirichlettrian)
  L = TrialFESpace(dirichlet_dofs[-vcat(fdofsd_new...)],M)

  # TO-DO: DIRTY (approach abandoned temporarily)
  # Begin compute DoF values Dirichlet via moment evaluation
  # (As far as I understand, only works for piece-wise constant functions on the facets,
  #  for higher degrees we need change of basis, or L2 projection)
  # dcvΓd_num=∫(mhscal*p)*dΓd
  # dcvΓd_den=∫(x->1.0)*dΓd
  # data_vΓd_num=Gridap.CellData.get_contribution(dcvΓd_num,dΓd.quad.trian)
  # data_vΓd_den=Gridap.CellData.get_contribution(dcvΓd_den,dΓd.quad.trian)
  # data_vΓd_num_den=lazy_map(Broadcasting(/),data_vΓd_num,data_vΓd_den)
  # dc = Gridap.CellData.DomainContribution()
  # Gridap.CellData.add_contribution!(dc, dirichlettrian, data_vΓd_num_den)
  # data = Gridap.FESpaces.collect_cell_vector(M,dc)
  # data[2] .= -data[2]
  # assem = SparseMatrixAssembler(M,L)
  # dirichlet_dof_values = zeros(num_dirichlet_dofs(M))
  # Gridap.FESpaces.assemble_vector!(dirichlet_dof_values,assem,data)
  # L=TrialFESpace(dirichlet_dof_values,M)

  # Define test FESpaces
  V = TestFESpace(model  , reffeᵤ; conformity=:L2)
  Q = TestFESpace(model  , reffeₚ; conformity=:L2)
  Y = MultiFieldFESpace([V,Q,M])

  # Create trial spaces
  U = TrialFESpace(V)
  P = TrialFESpace(Q)
  X = MultiFieldFESpace([U, P, L])

  yh = get_fe_basis(Y)
  vh,qh,mh = yh

  xh = get_trial_fe_basis(X)
  uh,ph,lh = xh

  trian = Triangulation(model)
  degree = 2*(order+1)
  dΩ = Measure(trian,degree)

  # neumanntags  = [5,6]
  # # TO-DO: neumanntrian = Triangulation(model_Γ,tags=neumanntags) this causes
  # # dcvΓ=∫(mh*g)*dΓ to fail in change_domain ...
  # neumanntrian = BoundaryTriangulation(model,tags=neumanntags)
  # degree = 2*(order+1)
  # dΓn = Measure(neumanntrian,degree)

  dcmΩ=∫( vh⋅uh - (∇⋅vh)*ph + qh*(∇⋅uh) )*dΩ
  dcvΩ=∫( vh⋅f + qh*(∇⋅u))*dΩ
  data_mΩ=Gridap.CellData.get_contribution(dcmΩ,dΩ.quad.trian)
  data_vΩ=Gridap.CellData.get_contribution(dcvΩ,dΩ.quad.trian)

  ∂T     = CellBoundary(model)
  x,w    = quadrature_evaluation_points_and_weights(∂T,2)

  #∫( mh*(uh⋅n) )*d∂K
  print("restrict_to_cell_boundary(∂T,uh)")
  @time uh_∂T = restrict_to_cell_boundary(∂T,uh)
  print("restrict_to_cell_boundary(∂T,mh)")
  @time mh_∂T = restrict_to_cell_boundary(∂T,mh)
  print("integrate_mh_mult_uh_cdot_n_low_level(∂T,mh_∂T,uh_∂T,x,w)")
  @time mh_mult_uh_cdot_n=integrate_mh_mult_uh_cdot_n_low_level(∂T,mh_∂T,uh_∂T,x,w)

  #∫( (vh⋅n)*lh )*d∂K
  print("restrict_to_cell_boundary(∂T,vh)")
  @time vh_∂T = restrict_to_cell_boundary(∂T,vh)
  print("restrict_to_cell_boundary(∂T,lh)")
  @time lh_∂T = restrict_to_cell_boundary(∂T,lh)
  print("integrate_vh_cdot_n_mult_lh_low_level(∂T,vh_∂T,lh_∂T,x,w)")
  @time vh_cdot_n_mult_lh=integrate_vh_cdot_n_mult_lh_low_level(∂T,vh_∂T,lh_∂T,x,w)

  print("Broadcasting(+),vh_cdot_n_mult_lh,mh_mult_uh_cdot_n,data_mΩ")
  @time cmat=lazy_map(Broadcasting(+),
                lazy_map(Broadcasting(+),vh_cdot_n_mult_lh,mh_mult_uh_cdot_n),
              data_mΩ)

  cvec=data_vΩ

  # cell=1
  # A11=vcat(hcat(cmat[cell][1,1],cmat[cell][1,2]),hcat(cmat[cell][2,1],0.0))
  # A12=vcat(cmat[cell][1,3],zeros(1,4))
  # A21=hcat(cmat[cell][3,1],zeros(4))
  # S22=-A21*inv(A11)*A12
  # b1=vcat(cvec[cell][1],cvec[cell][2])
  # y2=-A21*inv(A11)*b1

  k=StaticCondensationMap([1,2],[3])
  print("lazy_map(StaticCondensationMap,cmat,cvec)")
  @time cmat_cvec_condensed=lazy_map(k,cmat,cvec)

  #fdofsn=get_cell_dof_ids(M,neumanntrian)
  #fdofsd=get_cell_dof_ids(M,dirichlettrian)

  fdofscb=restrict_facet_dof_ids_to_cell_boundary(∂T,get_cell_dof_ids(M))
  assem = SparseMatrixAssembler(M,L)

  # This array is stored as a member variable of UnconstrainedFESpace in Gridap
  cell_is_dirichlet=Gridap.FESpaces._generate_cell_is_dirichlet(fdofscb)

  # This would be computed using get_cell_dof_values(lhd,cell_boundary) inside Gridap
  posnegreindexk=Gridap.Fields.PosNegReindex(zeros(num_free_dofs(M)),L.dirichlet_values)
  lhₖ= lazy_map(Gridap.Fields.Broadcasting(posnegreindexk),fdofscb)

  # This function is called from function _pair_contribution_when_possible(biform,liform,uhd)
  print("attach_dirichlet(cmat_cvec_condensed,lhₖ,cell_is_dirichlet)")
  @time cmat_cvec_condensed_dirichlet=Gridap.FESpaces.attach_dirichlet(cmat_cvec_condensed,
                                                                       lhₖ,
                                                                       cell_is_dirichlet)

  print("assemble_matrix_and_vector")
  @time A,b=assemble_matrix_and_vector(assem,
        (([cmat_cvec_condensed_dirichlet], [fdofscb], [fdofscb]),
                                            ([],[],[]),
                                            ([],[])))
#@time A,b=assemble_matrix_and_vector(assem,(([cmat_cvec_condensed], [fdofscb], [fdofscb]),
#                                      ([data_mΓd],[fdofsd],[fdofsd]),
#                                      ([data_vΓn,data_vΓd],[fdofsn,fdofsd])))
#@time A,b=assemble_matrix_and_vector(assem,(([cmat_cvec_condensed], [fdofscb], [fdofscb]),
#                                      ([data_mΓd],[fdofsd],[fdofsd]),
#                                      ([data_vΓd],[fdofsd])))
  x     = A\b
  lh    = FEFunction(L,x)

  k=BackwardStaticCondensationMap([1,2],[3])
  lhₖ= lazy_map(Gridap.Fields.Broadcasting(Gridap.Fields.PosNegReindex(
                         Gridap.FESpaces.get_free_dof_values(lh),lh.dirichlet_values)),
                      fdofscb)

  print("lazy_map(BackwardStaticCondensationMap,cmat,cvec,lhₖ)")
  @time uhphlhₖ=lazy_map(k,cmat,cvec,lhₖ)

  # tol=1.0e-12
  # cell=1
  # A11=vcat(hcat(cmat[cell][1,1],cmat[cell][1,2]),hcat(cmat[cell][2,1],0.0))
  # A12=vcat(cmat[cell][1,3],zeros(1,4))
  # A21=hcat(cmat[cell][3,1],zeros(4))

  # Am=vcat(hcat(A11,A12),hcat(A21,zeros(4,4)))
  # Am[6,6]=1.0
  # Am[7,7]=1.0
  # Am[8,8]=1.0
  # Am[9,9]=1.0
  # Sm=Am[6:9,6:9]-Am[6:9,1:5]*inv(Am[1:5,1:5])*Am[1:5,6:9]
  # ym=vcat(data_vΓd...)-A21*inv(Am[1:5,1:5])*vcat(cvec[1][1],cvec[1][2])
  # xm=Sm\ym
  # @assert norm(xm-x) < tol
  # bm=vcat(cvec[1][1],cvec[1][2],data_vΓd...)
  # Am\bm

  lhₑ=lazy_map(Gridap.Fields.BlockMap(3,3),ExploringGridapHybridization.convert_cell_wise_dofs_array_to_facet_dofs_array(∂T,
      lhₖ,get_cell_dof_ids(M)))

  assem = SparseMatrixAssembler(Y,X)
  lhₑ_dofs=get_cell_dof_ids(X,Triangulation(model_Γ))

  uhph_dofs=get_cell_dof_ids(X,Triangulation(model))
  uhph_dofs = lazy_map(Gridap.Fields.BlockMap(2,[1,2]),uhph_dofs.args[1],uhph_dofs.args[2])

  print("uh=lazy_map(x->x[1],uhphlhₖ)")
  @time uh=lazy_map(x->x[1],uhphlhₖ)
  print("ph=lazy_map(x->x[2],uhphlhₖ)")
  @time ph=lazy_map(x->x[2],uhphlhₖ)
  print("uhphₖ=lazy_map(Gridap.Fields.BlockMap(2,[1,2]),uh,ph)")
  @time uhphₖ=lazy_map(Gridap.Fields.BlockMap(2,[1,2]),uh,ph)

  print("free_dof_values=assemble_vector(assem,([lhₑ,uhphₖ],[lhₑ_dofs,uhph_dofs]))")
  @time free_dof_values=assemble_vector(assem,([lhₑ,uhphₖ],[lhₑ_dofs,uhph_dofs]))
  xh=FEFunction(X,free_dof_values)

  # cell=1
  # A11=vcat(hcat(cmat[cell][1,1],cmat[cell][1,2]),hcat(cmat[cell][2,1],0.0))
  # A12=vcat(cmat[cell][1,3],zeros(1,4))
  # A21=hcat(cmat[cell][3,1],zeros(4))
  # b1=vcat(cvec[cell][1],cvec[cell][2])
  # x2=lhₖ[cell]
  # x1=inv(A11)*(b1-A12*x2)
  #end
end

domain = (0,1,0,1)
partition = (1,2)
order = 0
model = CartesianDiscreteModel(domain,partition)
print("solve_darcy_rt_hdiv ")
@time sol_conforming=solve_darcy_rt_hdiv(model,order)
print("solve_darcy_hybrid_rt 1")
@time sol_nonconforming=solve_darcy_hybrid_rt(model,order)
print("solve_darcy_hybrid_rt 2")
@time sol_nonconforming=solve_darcy_hybrid_rt(model,order)
print("solve_darcy_hybrid_rt 3")
@time sol_nonconforming=solve_darcy_hybrid_rt(model,order)
trian = Triangulation(model)
degree = 2*(order+1)
dΩ = Measure(trian,degree)
uhc,_=sol_conforming
uhnc,_,_=sol_nonconforming

@test sqrt(sum(∫((uhc-uhnc)⋅(uhc-uhnc))dΩ)) < 1.0e-12

end # module
