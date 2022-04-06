
struct TransformOrthogonalDofBasis{Dc,Dp} <: Gridap.Arrays.Map end ;
struct OrthogonalBasisPullBack <: Gridap.Arrays.Map end;

function Gridap.FESpaces.get_cell_dof_basis(model::DiscreteModel,
       cell_reffe::AbstractArray{<:Gridap.ReferenceFEs.GenericRefFE{OrthogonalBasis}},
       ::Gridap.ReferenceFEs.L2Conformity)
    cell_map  = get_cell_map(Triangulation(model))
    phi       = cell_map[1]
    Jt        = lazy_map(Broadcasting(∇),cell_map)
    x         = lazy_map(get_nodes,lazy_map(get_dof_basis,cell_reffe))
    Jtx       = lazy_map(evaluate,Jt,x)
    reffe     = cell_reffe[1]
    Dc        = num_dims(reffe)
    et        = return_type(get_prebasis(reffe))
    pt        = Point{Dc,et}
    Dp        = first(size(return_type(phi,zero(pt))))
    k         = TransformOrthogonalDofBasis{Dc,Dp}()
    lazy_map(k,cell_reffe,Jtx)
end

function  Gridap.FESpaces.get_cell_shapefuns(model::DiscreteModel,
                            cell_reffe::AbstractArray{<:Gridap.ReferenceFEs.GenericRefFE{OrthogonalBasis}},
                            ::Gridap.ReferenceFEs.L2Conformity)
    cell_reffe_shapefuns=lazy_map(get_shapefuns,cell_reffe)
    k=OrthogonalBasisPullBack()
    lazy_map(k,
             cell_reffe_shapefuns,
             get_cell_map(Triangulation(model)))
end


function Gridap.Arrays.return_cache(::TransformOrthogonalDofBasis{Dc,Dp},
                                    reffe::Gridap.ReferenceFEs.GenericRefFE{OrthogonalBasis},
                                    Jtx) where {Dc,Dp}
  p = get_polytope(reffe)
  prebasis = get_prebasis(reffe)
  order = get_order(prebasis)
  et = return_type(prebasis)
  dofs = get_dof_basis(reffe)
  nodes, nf_nodes, nf_moments =  get_nodes(dofs),
                                 get_face_nodes_dofs(dofs),
                                 get_face_moments(dofs)
  db = MomentBasedDofBasis(nodes,nf_moments,nf_nodes)
  face_moments = [ similar(i)  for i in nf_moments ]

  cache = (db.nodes, db.face_nodes, nf_moments, face_moments)
  cache
end

function Gridap.Arrays.evaluate!(cache,
                                 ::TransformOrthogonalDofBasis,
                                 reffe::Gridap.ReferenceFEs.GenericRefFE{OrthogonalBasis},
                                 Jt_q)
  nodes, nf_nodes, nf_moments, face_moments = cache
  face_own_dofs=get_face_own_dofs(reffe)
  for face in 1:length(face_moments)
    nf_moments_face   = nf_moments[face]
    face_moments_face = face_moments[face]
    if length(nf_moments_face) > 0
      num_qpoints, num_moments = size(nf_moments_face)
      for i in 1:num_qpoints
        Jt_q_i = Jt_q[nf_nodes[face][i]]
        change = Gridap.ReferenceFEs.meas(Jt_q_i)
        for j in 1:num_moments
          face_moments_face[i,j] = change*nf_moments_face[i,j]
        end
      end
    end
  end
  MomentBasedDofBasis(nodes,face_moments,nf_nodes)
end

function Gridap.Arrays.evaluate!(cache,:: OrthogonalBasisPullBack,
                   v     :: Number,
                   detJ  :: Number)
  v*(1.0/detJ)
end

function Gridap.Arrays.evaluate!(cache,
                   k::OrthogonalBasisPullBack,
                   v::AbstractVector{<:Gridap.Fields.Field},
                   phi::Gridap.Fields.Field)
  Jt = ∇(phi)
  detJ = Operation(Gridap.ReferenceFEs.meas)(Jt)
  Broadcasting(Operation(k))(v,detJ)
end

function Gridap.Arrays.lazy_map(
  k::OrthogonalBasisPullBack,
  cell_ref_shapefuns::AbstractArray{<:AbstractArray{<:Gridap.Fields.Field}},
  cell_map::AbstractArray{<:Gridap.Fields.Field})

  cell_Jt = lazy_map(∇,cell_map)
  cell_detJ = lazy_map(Operation(Gridap.ReferenceFEs.meas),cell_Jt)

  lazy_map(Broadcasting(Operation(k)),cell_ref_shapefuns,cell_detJ)
end

function Gridap.Arrays.evaluate!(
  cache,
  ::Broadcasting{typeof(∇)},
  a::Gridap.Fields.BroadcastOpFieldArray{OrthogonalBasisPullBack})
  v, Jt = a.args
  # Assuming J comes from an affine map
  ∇v = Broadcasting(∇)(v)
  k = OrthogonalBasisPullBack()
  Broadcasting(Operation(k))(∇v,Jt)
end

function Gridap.Arrays.lazy_map(
  ::Broadcasting{typeof(gradient)},
  a::LazyArray{<:Fill{Broadcasting{Operation{OrthogonalBasisPullBack}}}})
  v, Jt = a.args
  ∇v = lazy_map(Broadcasting(∇),v)
  k = OrthogonalBasisPullBack()
  lazy_map(Broadcasting(Operation(k)),∇v,Jt)
end
