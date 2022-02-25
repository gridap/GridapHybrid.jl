
# These codes should go eventually to Gridap at the end!!!


function Gridap.CellData.get_triangulation(f::Gridap.MultiField.MultiFieldCellField)
  s1 = first(f.single_fields)
  trian = get_triangulation(s1)
  trian
end

##

function Gridap.CellData.CellQuadrature(trian::Triangulation,cell_quad,ids::DomainStyle)
  ctype_to_quad, cell_to_ctype = compress_cell_data(cell_quad)
  ctype_to_point = map(get_coordinates,ctype_to_quad)
  ctype_to_weight = map(get_weights,ctype_to_quad)
  cell_point = expand_cell_data(ctype_to_point,cell_to_ctype)
  cell_weight = expand_cell_data(ctype_to_weight,cell_to_ctype)
  CellQuadrature(cell_quad,cell_point,cell_weight,trian,ReferenceDomain(),ids)
end

function _get_f(v::Gridap.Fields.VectorBlock,f::Function)
  n  = length(v.array)
  v1 = f(v.array[1])
  TD = typeof(v1)
  a  = Vector{TD}(undef,n)
  t  = v.touched
  for i=1:n
    a[i]=f(v.array[i])
  end
  Gridap.Fields.ArrayBlock(a,t)
end

function Gridap.Geometry.get_coordinates(v::Gridap.Fields.VectorBlock)
  _get_f(v,Gridap.Geometry.get_coordinates)
end

function Gridap.Geometry.get_weights(v::Gridap.Fields.VectorBlock)
  _get_f(v,Gridap.Geometry.get_weights)
end

## Polynomials

# Optimizing evaluation at a single point

function Gridap.Arrays.return_cache(f::Gridap.Polynomials.QCurlGradMonomialBasis{D,T},x::Point) where {D,T}
  ndof = size(f.qgrad)[1]
  r = Gridap.Arrays.CachedArray(zeros(VectorValue{D,T},(ndof,)))
  xs = [x]
  cf = Gridap.Arrays.return_cache(f,xs)
  r, cf, xs
end

function Gridap.Arrays.evaluate!(cache,f::Gridap.Polynomials.QCurlGradMonomialBasis{D,T},x::Point) where {D,T}
  r, cf, xs = cache
  xs[1] = x
  v = Gridap.Arrays.evaluate!(cf,f,xs)
  ndof = size(v,2)
  Gridap.Arrays.setsize!(r,(ndof,))
  a = r.array
  copyto!(a,v)
  a
end

function Gridap.Arrays.return_cache(
  k::Gridap.Fields.BroadcastingFieldOpMap,
  f::Gridap.Fields.ArrayBlock{A,1},
  g::Gridap.Fields.ArrayBlock{B,1}) where {A,B}
  Gridap.Helpers.@check size(f.array) == (1,) ||
                        size(g.array) == (1,) ||
                        size(f.array) == size(g.array)
  fi = testvalue(A)
  gi = testvalue(B)
  ci = return_cache(k,fi,gi)
  hi = evaluate!(ci,k,fi,gi)
  s = (max(size(f.array,1),size(g.array,1)),)
  a = Array{typeof(hi),1}(undef,s)
  b = Array{typeof(ci),1}(undef,s)

  m = Gridap.Fields.ZeroBlockMap()
  zf = Array{typeof(return_cache(m,fi,gi))}(undef,size(f.array))
  zg = Array{typeof(return_cache(m,gi,fi))}(undef,size(f.array))

  t = fill(false,s)
  for i=1:s[1]
     findex=min(i,size(f.array)[1])
     gindex=min(i,size(g.array)[1])
     if f.touched[findex] && g.touched[gindex]
       b[i] = return_cache(k,f.array[findex],g.array[gindex])
       t[i] = true
     else
        if (size(f.array)==size(g.array))
           if f.touched[findex]
            _fi = f.array[findex]
            zg[i] = return_cache(m,gi,_fi)
            _gi = evaluate!(zg[i],m,gi,_fi)
            b[i] = return_cache(k,_fi,_gi)
            t[i] = true
           elseif g.touched[gindex]
            _gi = g.array[gindex]
            zf[i] = return_cache(m,fi,_gi)
            _fi = evaluate!(zf[i],m,fi,_gi)
            b[i] = return_cache(k,_fi,_gi)
            t[i] = true
           end
        end
     end
  end
  ArrayBlock(a,t), b, zf, zg
end

function Gridap.Arrays.evaluate!(
  cache,
  k::Gridap.Fields.BroadcastingFieldOpMap,
  f::Gridap.Fields.ArrayBlock{A,1},
  g::Gridap.Fields.ArrayBlock{B,1}) where {A,B}
  Gridap.Helpers.@check size(f.array) == (1,) ||
                        size(g.array) == (1,) ||
                        size(f.array) == size(g.array)
  a,b,zf,zg=cache
  s = (max(size(f.array,1),size(g.array,1)),)
  m = Gridap.Fields.ZeroBlockMap()
  for i=1:s[1]
     findex=min(i,size(f.array)[1])
     gindex=min(i,size(g.array)[1])
     if f.touched[findex] && g.touched[gindex]
        Gridap.Helpers.@check a.touched[i]
        a.array[i] = evaluate!(b[i],k,f.array[findex],g.array[gindex])
     else
        if (size(f.array)==size(g.array))
          if f.touched[findex]
              Gridap.Helpers.@check a.touched[i]
              fi = f.array[findex]
              gi = evaluate!(zg[i],m,nothing,fi)
              a.array[i] = evaluate!(b[i],k,fi,gi)
          elseif g.touched[gindex]
              Gridap.Helpers.@check a.touched[i]
              gi = g.array[gindex]
              fi = evaluate!(zf[i],m,nothing,gi)
              a.array[i] = evaluate!(b[i],k,fi,gi)
          end
        end
     end
  end
  a
end

function Gridap.Arrays.return_cache(
  k::Gridap.Fields.BroadcastingFieldOpMap,
  f::Gridap.Fields.ArrayBlock{A,1},
  g::Gridap.Fields.ArrayBlock{B,2}) where {A,B}
  # Degenerated case in which we have a single-block f
  if (size(f)==(1,))
    return return_cache(k,f.array[1],g)
  end
  fi = testvalue(A)
  gi = testvalue(B)
  ci = return_cache(k,fi,gi)
  hi = evaluate!(ci,k,fi,gi)
  Gridap.Helpers.@check size(g.array,1) == 1 || size(g.array,2) == 0
  s = (size(f.array,1),size(g.array,2))
  a = Array{typeof(hi),2}(undef,s)
  b = Array{typeof(ci),2}(undef,s)
  t = fill(false,s)
  for j in 1:s[2]
    for i in 1:s[1]
      if f.touched[i] && g.touched[1,j]
        t[i,j] = true
        b[i,j] = return_cache(k,f.array[i],g.array[1,j])
      end
    end
  end
  ArrayBlock(a,t), b
end

function Gridap.Arrays.evaluate!(
  cache,k::Gridap.Fields.BroadcastingFieldOpMap,
  f::Gridap.Fields.ArrayBlock{A,1},
  g::Gridap.Fields.ArrayBlock{B,2}) where {A,B}
  # Degenerated case in which we have a single-block f
  if (size(f)==(1,))
    return evaluate!(cache,k,f.array[1],g)
  end
  a,b = cache
  s = size(a.array)
  for j in 1:s[2]
    for i in 1:s[1]
      if f.touched[i] && g.touched[1,j]
        a.array[i,j] = evaluate!(b[i,j],k,f.array[i],g.array[1,j])
      end
    end
  end
  a
end

### The functions in the sequel are required to be able to compute jacobians of
### residuals composed of integrals over the skeleton of the mesh


function Gridap.FESpaces._jacobian(f,
  uh::Gridap.MultiField.MultiFieldFEFunction,
  fuh::Gridap.CellData.DomainContribution)
  terms = DomainContribution()
  U = Gridap.FESpaces.get_fe_space(uh)
  fuh = _merge_bulk_and_skeleton_vector_contributions(fuh)
  for trian in get_domains(fuh)
    Gridap.Helpers.@check isa(trian,SkeletonTriangulation)
    g = Gridap.FESpaces._change_argument(jacobian,f,trian,uh)
    cell_u = lazy_map(DensifyInnerMostBlockLevelMap(),get_cell_dof_values(uh,trian))
    cell_id = Gridap.FESpaces._compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_jacobian(g,cell_u,cell_id)
    monolithic_result=cell_grad
    blocks        = [] # TO-DO type unstable. How can I infer the type of its entries?
    blocks_coords = Tuple{Int,Int}[]
    nfields = length(U.spaces)
    cell_dofs_field_offsets=Gridap.MultiField._get_cell_dofs_field_offsets(uh,trian)
    for j=1:nfields
      view_range_j=cell_dofs_field_offsets[j]:cell_dofs_field_offsets[j+1]-1
      for i=1:nfields
        view_range_i=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
        # TO-DO: depending on the residual being differentiated, we may end with
        #        blocks [i,j] full of zeros. I guess that it might desirable to early detect
        #        these zero blocks and use a touch[i,j]==false block in ArrayBlock.
        #        How can we detect that we have a zero block?
        block=lazy_map(x->view(x,view_range_i,view_range_j),monolithic_result)
        append!(blocks,[block])
        append!(blocks_coords,[(i,j)])
      end
    end
    cell_grad=lazy_map(BlockMap((nfields,nfields),blocks_coords),blocks...)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end


# TO IMPROVE .... Perhaps via a Map?
# Do we have already this function in Gridap?
# What if there are no cells? I am assuming that there is at least one cell
# What if the number of dofs per field per cell is different among cells?
function Gridap.MultiField._get_cell_dofs_field_offsets(
  uh::Gridap.MultiField.MultiFieldFEFunction,trian::SkeletonTriangulation)
  U = Gridap.FESpaces.get_fe_space(uh)
  uh_dofs = get_cell_dof_values(uh,trian)[1]
  nfields = length(U.spaces)
  dofs_field_offsets=Vector{Int}(undef,nfields+1)
  dofs_field_offsets[1]=1
  for i in 1:nfields
    dofs_field_offsets[i+1]=dofs_field_offsets[i]+length(uh_dofs.array[i])
  end
  dofs_field_offsets
end

function _get_cell_dofs_facet_offsets(cell_lfacet_dofs::ArrayBlock)
  U = Gridap.FESpaces.get_fe_space(uh)
  uh_dofs = get_cell_dof_values(uh,trian)[1]
  nfields = length(U.spaces)
  dofs_field_offsets=Vector{Int}(undef,nfields+1)
  dofs_field_offsets[1]=1
  for i in 1:nfields
    dofs_field_offsets[i+1]=dofs_field_offsets[i]+length(uh_dofs.array[i])
  end
  dofs_field_offsets
end


struct SkeletonVectorFromNonBlockedSkeletonVector{A} <: Gridap.Arrays.Map
  cell_wise_block_sizes::A
end

function Gridap.Arrays.return_cache(k::SkeletonVectorFromNonBlockedSkeletonVector,
                                    v::AbstractVector,
                                    cell::Integer)
  #TO-DO: generate CompressedArray cache with all cell types
  cwbs_cache=array_cache(k.cell_wise_block_sizes)
  cwbsi=getindex!(cwbs_cache,k.cell_wise_block_sizes,cell)
  nb=length(cwbsi)
  touched=Vector{Bool}(undef,nb)
  touched.=true
  T=typeof(view(v,1:length(v)))
  array=Vector{T}(undef,nb)
  ArrayBlock(array,touched),cwbs_cache
end


function Gridap.Arrays.evaluate!(cache,
  k::SkeletonVectorFromNonBlockedSkeletonVector,
  v::AbstractVector,
  cell::Integer)
  output,cwbs_cache=cache
  cwbsi=getindex!(cwbs_cache,k.cell_wise_block_sizes,cell)
  nb=length(cwbsi)
  Gridap.Helpers.@check nb==length(output)
  current=1
  next=1
  for i=1:nb
    next=current+cwbsi[i]
    output.array[i]=view(v,current:next-1)
    current=next
  end
  output
end

function _generate_cell_lface_dofs_from_cell_dofs(glue,
                                                  cell_wise_facets,
                                                  cell_values_field,
                                                  facets_dofs_ids)
  facet_sizes=lazy_map(x->length(x),facets_dofs_ids)
  cell_wise_facet_sizes=SkeletonVectorFromFacetVector(glue,cell_wise_facets,facet_sizes)
  m=SkeletonVectorFromNonBlockedSkeletonVector(cell_wise_facet_sizes)
  lazy_map(m,cell_values_field,collect(1:length(cell_values_field)))
end



function Gridap.FESpaces._change_argument(
  op::typeof(jacobian),f,trian::SkeletonTriangulation,uh::Gridap.MultiField.MultiFieldFEFunction)

  U = Gridap.FESpaces.get_fe_space(uh)
  function g(cell_u)
    single_fields = GenericCellField[]
    nfields = length(U.spaces)
    cell_dofs_field_offsets=Gridap.MultiField._get_cell_dofs_field_offsets(uh,trian)
    for i in 1:nfields
      view_range=cell_dofs_field_offsets[i]:cell_dofs_field_offsets[i+1]-1
      cell_values_field = lazy_map(a->view(a,view_range),cell_u)
      Ui=U.spaces[i]
      Ui_trian=get_triangulation(Ui)
      model=get_background_model(Ui_trian)
      cell_wise_facets = _get_cell_wise_facets(model)
      D=num_cell_dims(model)
      if (isa(Ui_trian, Triangulation{D,D}))         # Cell triangulation
        cf = CellField(Ui, cell_values_field)
      elseif (isa(Ui_trian, Triangulation{D - 1,D})) # Facet triangulation
        cell_lface_dof_values =
         _generate_cell_lface_dofs_from_cell_dofs(
            trian.glue,
            cell_wise_facets,
            cell_values_field,
            get_cell_dof_ids(Ui))
        Ui_basis = get_fe_basis(Ui)
        Ui_basis_data = get_data(Ui_basis)
        Ui_basis_cell_lface_data = SkeletonVectorFromFacetVector(
          trian.glue,
          cell_wise_facets,
          Ui_basis_data)
        cell_field = lazy_map(linear_combination, cell_lface_dof_values, Ui_basis_cell_lface_data)
        cf = GenericCellField(cell_field, trian, DomainStyle(Ui_basis))
      else
        Gridap.Helpers.@notimplemented
      end
      push!(single_fields,cf)
    end
    xh = Gridap.MultiField.MultiFieldCellField(single_fields)
    cell_grad = f(xh)
    cell_grad=_merge_bulk_and_skeleton_matrix_contributions(cell_grad)
    cell_grad_cont_block=get_contribution(cell_grad,trian)
    bs = [cell_dofs_field_offsets[i+1]-cell_dofs_field_offsets[i] for i=1:nfields]
    lazy_map(DensifyInnerMostBlockLevelMap(),
             Fill(bs,length(cell_grad_cont_block)),
             cell_grad_cont_block)
  end
  g
end
