module DarcyRTHTests

using Test
using Gridap
using FillArrays
using Gridap.Geometry
using ExploringGridapHybridization

#2D problem
u2(x) = VectorValue(1+x[1],1+x[2])
Gridap.divergence(::typeof(u2)) = (x) -> 2
p(x) = -3.14
∇p(x) = zero(x)
Gridap.∇(::typeof(p)) = ∇p
f2(x) = u2(x) + ∇p(x)
# Normal component of u(x) on Neumann boundary
function g2(x)
  tol=1.0e-14
  if (abs(x[2])<tol)
    return -x[2]
  elseif (abs(x[2]-1.0)<tol)
    return x[2]
  end
  Gridap.Helpers.@check false
end

#3D problem
u3(x) = VectorValue(1+x[1],1+x[2],1+x[3])
Gridap.divergence(::typeof(u3)) = (x) -> 3
f3(x) = u3(x) + ∇p(x)
function g3(x) # Normal component of u(x) on Neumann boundary
  @assert false
end

function ufg(D::Int)
   if (D==2)
    u2,f2,g2
   elseif (D==3)
    u3,f3,g3
   end
end

function dirichlet_tags(D::Int)
  if (D==2)
    collect(5:8)
  elseif (D==3)
    collect(21:26)
  end
end

function pre_assembly_stage_rt_hdiv(model,order)
  D=num_cell_dims(model)
  dtags=dirichlet_tags(D)
  u,f,_ = ufg(D)
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
  btrian = BoundaryTriangulation(model,tags=dtags)
  dΓ = Measure(btrian,degree)
  nb = get_normal_vector(btrian)
  a((u, p),(v, q)) = ∫( u⋅v - (∇⋅v)*p + q*(∇⋅u) )*dΩ
  b(( v, q)) = ∫( v⋅f + q*(∇⋅u))*dΩ - ∫((v⋅nb)*p )*dΓ
  X,Y,a,b
end

function assembly_stage_rt_hdiv(X,Y,a,b)
  AffineFEOperator(a,b,X,Y)
end

function solve_stage_rt_hdiv(op)
  xh = solve(op)
end

function solve_darcy_rt_hdiv(model,order)
  X,Y,a,b=pre_assembly_stage_rt_hdiv(model,order)
  op=assembly_stage_rt_hdiv(X,Y,a,b)
  solve_stage_rt_hdiv(op)
end

macro mytime(expr,concept,activate)
  if (activate)
   print(concept)
   @time esc(expr)
  else
   esc(expr)
  end
end

function preassembly_stage_darcy_hybrid_rt(model,∂T,order)
  # Geometry part
  D=num_cell_dims(model)
  model_Γ = BoundaryDiscreteModel(Polytope{D-1},model,collect(1:num_facets(model)))
  u,f,_ = ufg(D)

  # Functional part
  # To investigate what is needed to have an inf-sup stable triplet
  # for the RT-H method
  reffeᵤ = ReferenceFE(raviart_thomas,Float64,order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order)
  reffeₗ = ReferenceFE(lagrangian,Float64,order)

  # Compute the Dof values of Dirichlet DoFs (L2 projection)
  M = TestFESpace(model_Γ,
                  reffeₗ;
                  conformity=:L2,
                  dirichlet_tags=num_faces(get_polytopes(model)[1]))
  L = TrialFESpace(M)
  dtags=dirichlet_tags(D)
  dirichlettrian=BoundaryTriangulation(model,tags=dtags)
  degree = 2*(order+1)
  dΓd = Measure(dirichlettrian,degree)
  mh = get_fe_basis(M)
  lh = get_trial_fe_basis(L)
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
  M = TestFESpace(model_Γ,
                  reffeₗ;
                  conformity=:L2,
                  dirichlet_tags=dtags)
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

  x,w    = quadrature_points_and_weights(∂T,2)

  #∫( mh*(uh⋅n) )*d∂K
  @mytime uh_∂T=restrict_to_cell_boundary(∂T,uh) "restrict_to_cell_boundary(∂T,uh)" false
  @mytime mh_∂T=restrict_to_cell_boundary(∂T,mh) "restrict_to_cell_boundary(∂T,mh)" false
  @mytime mh_mult_uh_cdot_n=integrate_mh_mult_uh_cdot_n_low_level(∂T,mh_∂T,uh_∂T,x,w) "integrate_mh_mult_uh_cdot_n_low_level(∂T,mh_∂T,uh_∂T,x,w)" false

  #∫( (vh⋅n)*lh )*d∂K
  @mytime vh_∂T = restrict_to_cell_boundary(∂T,vh) "restrict_to_cell_boundary(∂T,vh)" false
  @mytime lh_∂T = restrict_to_cell_boundary(∂T,lh) "restrict_to_cell_boundary(∂T,lh)" false
  @mytime vh_cdot_n_mult_lh=integrate_vh_cdot_n_mult_lh_low_level(∂T,vh_∂T,lh_∂T,x,w) "integrate_vh_cdot_n_mult_lh_low_level(∂T,vh_∂T,lh_∂T,x,w)" false

  @mytime cmat=lazy_map(Broadcasting(+),lazy_map(Broadcasting(+),vh_cdot_n_mult_lh,mh_mult_uh_cdot_n),data_mΩ) "Broadcasting(+),vh_cdot_n_mult_lh,mh_mult_uh_cdot_n,data_mΩ" false

  cvec=data_vΩ
  k=StaticCondensationMap([1,2],[3])
  @mytime cmat_cvec_condensed=lazy_map(k,cmat,cvec) "lazy_map(StaticCondensationMap,cmat,cvec)" false
  model_Γ,X,Y,M,L,cmat,cvec,cmat_cvec_condensed
end

function assembly_stage_darcy_hybrid_rt(model,∂T,M,L,cmat_cvec_condensed)
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
  @mytime cmat_cvec_condensed_dirichlet=Gridap.FESpaces.attach_dirichlet(cmat_cvec_condensed,lhₖ,cell_is_dirichlet) "attach_dirichlet(cmat_cvec_condensed,lhₖ,cell_is_dirichlet)" false

  @mytime A,b=assemble_matrix_and_vector(assem,(([cmat_cvec_condensed_dirichlet], [fdofscb], [fdofscb]),([],[],[]),([],[]))) "assemble_matrix_and_vector" false
  #@mytime A,b=assemble_matrix_and_vector(assem,(([cmat_cvec_condensed], [fdofscb], [fdofscb]),
  #                                      ([data_mΓd],[fdofsd],[fdofsd]),
  #                                      ([data_vΓn,data_vΓd],[fdofsn,fdofsd])))
  #@mytime A,b=assemble_matrix_and_vector(assem,(([cmat_cvec_condensed], [fdofscb], [fdofscb]),
  #                                      ([data_mΓd],[fdofsd],[fdofsd]),
  #                                      ([data_vΓd],[fdofsd])))
  # cell=1
  # A11=vcat(hcat(cmat[cell][1,1],cmat[cell][1,2]),hcat(cmat[cell][2,1],0.0))
  # A12=vcat(cmat[cell][1,3],zeros(1,4))
  # A21=hcat(cmat[cell][3,1],zeros(4))
  # S22=-A21*inv(A11)*A12
  # b1=vcat(cvec[cell][1],cvec[cell][2])
  # y2=-A21*inv(A11)*b1
  A,b,fdofscb
end

function solve_stage_darcy_hybrid_rt(A,b)
  x = A\b
end

function back_substitution_stage_darcy_hybrid_rt(∂T,model,model_Γ,X,Y,M,L,x,cmat,cvec,fdofscb)
  lh    = FEFunction(L,x)

  k=BackwardStaticCondensationMap([1,2],[3])
  lhₖ= lazy_map(Gridap.Fields.Broadcasting(Gridap.Fields.PosNegReindex(
                         Gridap.FESpaces.get_free_dof_values(lh),lh.dirichlet_values)),
                      fdofscb)
  @mytime uhphlhₖ=lazy_map(k,cmat,cvec,lhₖ) "lazy_map(BackwardStaticCondensationMap,cmat,cvec,lhₖ)" false

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

  @mytime uhphₖ=lazy_map(RestrictArrayBlockMap([1,2]),uhphlhₖ) "uhphₖ=lazy_map(RestrictArrayBlockMap([1,2]),uhphlhₖ)" false

  @mytime free_dof_values=assemble_vector(assem,([lhₑ,uhphₖ],[lhₑ_dofs,uhph_dofs])) "free_dof_values=assemble_vector(assem,([lhₑ,uhphₖ],[lhₑ_dofs,uhph_dofs]))" false
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

function solve_darcy_hybrid_rt(model,∂T,order)
  model_Γ,X,Y,M,L,cmat,cvec,cmat_cvec_condensed=preassembly_stage_darcy_hybrid_rt(model,∂T,order)
  A,b,fdofscb=assembly_stage_darcy_hybrid_rt(model,∂T,M,L,cmat_cvec_condensed)
  x=solve_stage_darcy_hybrid_rt(A,b)
  back_substitution_stage_darcy_hybrid_rt(∂T,model,model_Γ,X,Y,M,L,x,cmat,cvec,fdofscb)
end

domain = (0,1,0,1)
partition = (2,2)
order = 0
model = CartesianDiscreteModel(domain,partition)

sol_conforming=solve_darcy_rt_hdiv(model,order)
∂T = CellBoundary(model)
sol_nonconforming=solve_darcy_hybrid_rt(model,∂T,order)
trian = Triangulation(model)
degree = 2*(order+1)
dΩ = Measure(trian,degree)
uhc,_=sol_conforming
uhnc,_,_=sol_nonconforming
@test sqrt(sum(∫((uhc-uhnc)⋅(uhc-uhnc))dΩ)) < 1.0e-12

end # module
