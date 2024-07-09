module OrthogonalBasisRefFEsTests
  using Gridap
  using GridapHybrid
  using Test
  using LinearAlgebra

  function orthogonal_basis_tests(model, T, order, subspace_trial, subspace_test)
    D          = num_cell_dims(model)
    reffe_trial = ReferenceFE(orthogonal_basis, T, order; subspace=subspace_trial)
    reffe_test  = ReferenceFE(monomial_basis, T, order; subspace=subspace_test)
    # We need to explicitly pass vector_type here so that we avoid that
    # TestFESpace constructor invokes UndefinedDofBasis
    V          = TestFESpace(model, reffe_test; conformity = :L2, vector_type=Vector{Float64})
    U          = TrialFESpace(TestFESpace(model, reffe_trial; conformity = :L2))
    Ω          = Triangulation(ReferenceFE{D}, model)
    degree     = 2 * (order + 1) # TO-DO: To think which is the minimum degree required
    dΩ         = Measure(Ω, degree)
    a(u,v)=∫(v⋅u)dΩ
    A=assemble_matrix(a,U,V)
    @test A ≈ I

    if (T==Float64 && subspace_trial==:ZeroMean)
      du=get_fe_basis(U)
      @test all( norm.(get_array(∫(du)dΩ)) .< 1.0e-14 )
    elseif (T==Float64 && subspace_trial==:NonZeroMean)
      du=get_fe_basis(U)
      @test all( norm.(get_array(∫(du)dΩ)) .≈ 1.0 )
    end

  end
  partition = (0,0.25,0,0.25)
  cells     = (10,10)
  model     = CartesianDiscreteModel(partition, cells)

  for T in (Float64,VectorValue{2,Float64})
     orthogonal_basis_tests(model,T,1,:Full,:Full)
     orthogonal_basis_tests(model,T,1,:ZeroMean,:ExcludeConstant)
     orthogonal_basis_tests(model,T,1,:NonZeroMean,:OnlyConstant)
  end

end
