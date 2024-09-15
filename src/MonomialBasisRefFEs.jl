struct UndefinedDofBasis <: AbstractVector{Gridap.ReferenceFEs.Dof} end;

Base.size(a::UndefinedDofBasis) = @assert false
Base.axes(a::UndefinedDofBasis) = @assert false
Base.getindex(a::UndefinedDofBasis,i::Integer) = @assert false
Base.IndexStyle(::UndefinedDofBasis) = @assert false

function Gridap.Arrays.return_cache(b::UndefinedDofBasis,field)
  @assert false
end

function Gridap.Arrays.evaluate!(cache,b::UndefinedDofBasis,field)
  @assert false
end



struct MonomialBasis <: Gridap.ReferenceFEs.ReferenceFEName end

const monomial_basis = MonomialBasis()

const monomial_basis_reffe_subspaces = (:Full,:OnlyConstant,:ExcludeConstant)

"""
    MonomialBasisRefFE(::Type{T},p::Polytope,order::Integer; subspace) where T
    Possible values for subspace are:
         * :Full           : Includes all monomials.
         * :OnlyConstant   : Subspace .
         * :ExcludeConstant: Subspace .
"""

_mb_p_filter_full(e,order)             = (sum(e) <= order)
_mb_p_filter_only_constant(e,order)    = (sum(e) == 0)
_mb_p_filter_exclude_constant(e,order) = (0 < sum(e) <= order)

_mb_p_filter = Dict(:Full            =>_mb_p_filter_full,
                    :OnlyConstant    =>_mb_p_filter_only_constant,
                    :ExcludeConstant =>_mb_p_filter_exclude_constant)

function MonomialBasisRefFE(::Type{T},p::Polytope,order::Integer;subspace=:Full) where T

  Gridap.Helpers.@check (subspace in monomial_basis_reffe_subspaces) """\n
  Possible values of subspace are $(monomial_basis_reffe_subspaces).
  """

  prebasis = Gridap.ReferenceFEs.MonomialBasis{num_dims(p)}(T,order,_mb_p_filter[subspace])

  face_dofs = [Int[] for i=1:num_faces(p)-1]
  push!(face_dofs, [i for i=1:length(prebasis)])
  dof_basis = UndefinedDofBasis()
  ndofs = length(prebasis)
  metadata = nothing
  conf = Gridap.ReferenceFEs.L2Conformity()

  reffe = Gridap.ReferenceFEs.GenericRefFE{MonomialBasis}(
    ndofs,
    p,
    prebasis,
    dof_basis,
    conf,
    metadata,
    face_dofs,
    prebasis)

  reffe
end

function Gridap.ReferenceFEs.ReferenceFE(p::Polytope,
                                          ::MonomialBasis,
                                          order;
                                          subspace=:Full)
  MonomialBasisRefFE(Float64,p,order;subspace=subspace)
end

function Gridap.ReferenceFEs.ReferenceFE(p::Polytope,
                                          ::MonomialBasis,
                                          ::Type{T},
                                          order;
                                          subspace=:Full) where T
  MonomialBasisRefFE(T,p,order;subspace=subspace)
end

function Gridap.ReferenceFEs.Conformity(reffe::Gridap.ReferenceFEs.GenericRefFE{MonomialBasis},sym::Symbol)
  if sym == :L2
    L2Conformity()
  else
    Gridap.Helpers.@unreachable """\n
    It is not possible to use conformity = $sym on a OrthogonalBasis reference FE.

    Possible values of conformity for this reference fe are $((:L2)).
    """
  end
end
