
struct OrthogonalBasis <: Gridap.ReferenceFEs.ReferenceFEName end

const orthogonal_basis = OrthogonalBasis()

const orthogonal_basis_reffe_subspaces = (:Full,:NonZeroMean,:ZeroMean)

"""
    OrthogonalBasisRefFE(::Type{T},p::Polytope,order::Integer; subspace) where T
    Possible values for subspace are:
         * :Full       : Includes all polynomials.
         * :NonZeroMean: Subspace of orthogonal polynomials with non-zero mean over the cell.
         * :ZeroMean   : Subspace of orthogonal polynomials with zero mean over the cell.
"""


function OrthogonalBasisRefFE(::Type{T},p::Polytope,order::Integer; subspace::Symbol=:Full) where T

  Gridap.Helpers.@check (subspace in orthogonal_basis_reffe_subspaces) """\n
  Possible values of subspace are $(orthogonal_basis_reffe_subspaces).
  """

  D = num_dims(p)
  extrusion = Gridap.ReferenceFEs.tfill(Gridap.ReferenceFEs.TET_AXIS,Val{D}())
  simplex = Gridap.ReferenceFEs.ExtrusionPolytope(extrusion)
  reffe = Gridap.ReferenceFEs.LagrangianRefFE(T,simplex,order)

  prebasis = Gridap.ReferenceFEs.get_prebasis(reffe)

  nodes, moments = _orthogonal_basis_nodes_and_moments(T,p,order)

  face_dofs = Gridap.ReferenceFEs._face_own_dofs_from_moments(moments)

  dof_basis = Gridap.ReferenceFEs.MomentBasedDofBasis(nodes, moments)

  shapefuns=Gridap.ReferenceFEs.compute_shapefuns(dof_basis,prebasis)
  shapefuns=_orthogonal_basis_filter_shapefuns(T,shapefuns,subspace)

  moments = _orthogonal_basis_filter_moments(T,moments,subspace)
  dof_basis = Gridap.ReferenceFEs.MomentBasedDofBasis(nodes, moments)

  ndofs = Gridap.ReferenceFEs.num_dofs(dof_basis)

  metadata = nothing

  conf = Gridap.ReferenceFEs.L2Conformity()

  reffe = Gridap.ReferenceFEs.GenericRefFE{OrthogonalBasis}(
    ndofs,
    p,
    prebasis,
    dof_basis,
    conf,
    metadata,
    face_dofs,
    shapefuns)

  reffe
end

function Gridap.ReferenceFEs.ReferenceFE(p::Polytope,
                                          ::OrthogonalBasis,
                                          order;
                                          subspace=:Full)
  OrthogonalBasisRefFE(Float64,p,order;subspace=subspace)
end

function Gridap.ReferenceFEs.ReferenceFE(p::Polytope,
                                          ::OrthogonalBasis,
                                          ::Type{T},
                                          order;
                                          subspace=:Full) where T
  OrthogonalBasisRefFE(T,p,order;subspace=subspace)
end

function Gridap.ReferenceFEs.Conformity(reffe::Gridap.ReferenceFEs.GenericRefFE{OrthogonalBasis},sym::Symbol)
  if sym == :L2
    L2Conformity()
  else
    Gridap.Helpers.@unreachable """\n
    It is not possible to use conformity = $sym on a OrthogonalBasis reference FE.

    Possible values of conformity for this reference fe are $((:L2)).
    """
  end
end

function _orthogonal_basis_moments(p, prebasis, ips, wips)
  # Interior DOFs-related basis evaluated at interior integration points
  ishfs_iips = evaluate(prebasis,ips)
  return wips.⋅ishfs_iips
end

function _orthogonal_basis_nodes_and_moments(::Type{T},
                                             p::Gridap.ReferenceFEs.Polytope,
                                             order::Integer) where T

  D = num_dims(p)
  ft = T
  et = eltype(T)
  pt = Point{D,et}

  nf_nodes = [ zeros(pt,0) for face in 1:num_faces(p)]
  nf_moments = [ zeros(ft,0,0) for face in 1:num_faces(p)]


  ccips, cmoments = _orthogonal_basis_nodes_and_moments(p,T,order)
  crange = Gridap.ReferenceFEs.get_dimrange(p,D)
  nf_nodes[crange] = ccips
  nf_moments[crange] = cmoments

  nf_nodes, nf_moments
end


_p_filter(e,order) = (sum(e) <= order)


# It provides for every cell the nodes and the moments arrays
function _orthogonal_basis_nodes_and_moments(p,T,order)
  degree = 2*order+1
  iquad  = Gridap.ReferenceFEs.Quadrature(p,degree)
  ips    = Gridap.ReferenceFEs.get_coordinates(iquad)
  wips   = Gridap.ReferenceFEs.get_weights(iquad)

  # Moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  prebasis = Gridap.ReferenceFEs.MonomialBasis{num_dims(p)}(T,order,_p_filter)
  moments = _orthogonal_basis_moments(p,prebasis,ips,wips)

  return [ips], [moments]

end


function _orthogonal_basis_filter_moments(T,moments,subspace)
  if subspace==:Full
     return moments
  else
     filtered_moments=[moments[i] for i=1:(length(moments)-1)]
     # All moments are associated to the cells in orthogonal basis RefFEs
     cmoments = moments[end]
     if subspace==:ZeroMean
      s = num_components(T)+1
      e = size(cmoments,2)
     elseif subspace==:NonZeroMean
      s = 1
      e = num_components(T)
     end
     push!(filtered_moments,cmoments[:,s:e])
     return filtered_moments
  end
end

function _orthogonal_basis_filter_shapefuns(T,
            shapefuns::Gridap.Fields.LinearCombinationFieldVector,subspace)
  if subspace==:Full
    return shapefuns
  else
    if subspace==:ZeroMean
      s = num_components(T)+1
      e = size(shapefuns.values,2)
     elseif subspace==:NonZeroMean
      s = 1
      e = num_components(T)
     end
     values=shapefuns.values[:,s:e]
     Gridap.Fields.LinearCombinationFieldVector(values,shapefuns.fields)
  end
end
