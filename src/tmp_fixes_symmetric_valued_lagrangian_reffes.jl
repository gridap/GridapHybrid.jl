
# These function overloads must go to Gridap 0.17, in the present (or improved) form !!!!

@generated function Gridap.TensorValues._flatten_upper_triangle(data::AbstractArray,::Val{D}) where D
  str = ""
  for i in 1:D
    for j in i:D
      str *= "data[$i,$j], "
    end
  end
  Meta.parse("($str)")
end

# length(T) is in turn used by num_components(T). I think that num_components(T) must be
# D*(D+1)÷2 and not D*D as it is defined now in Gridap v0.17.7
Base.length(::Type{<:Gridap.TensorValues.SymTensorValue{D}}) where {D} = D*(D+1)÷2

function Gridap.Polynomials._set_value!(
    v::AbstractVector{<:Gridap.TensorValues.SymTensorValue{D}},s::T,k) where {D,T}
    V = eltype(v)
    m = Vector{T}(undef,length(V)) # Allocation!
    z = zero(T)
    js = eachindex(m)
    for l in js
      for i in js
        @inbounds m[i] = z
      end
      m[l]=s
      v[k]=Tuple(m)
      k += 1
    end
    k
end

function Gridap.Polynomials._gradient_nd!(
  v::AbstractVector{G},
  x,
  orders,
  terms::AbstractVector{CartesianIndex{D}},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{<:Gridap.TensorValues.SymTensorValue{D}}) where {G,T,D}

  dim = D
  for d in 1:dim
    Gridap.Polynomials._evaluate_1d!(c,x,orders[d],d)
    Gridap.Polynomials._gradient_1d!(g,x,orders[d],d)
  end

  z = zero(Gridap.TensorValues.Mutable(VectorValue{D,T}))
  o = one(T)
  k = 1

  for ci in terms

    s = z
    for i in eachindex(s)
      @inbounds s[i] = o
    end
    for q in 1:dim
      for d in 1:dim
        if d != q
          @inbounds s[q] *= c[d,ci[d]]
        else
          @inbounds s[q] *= g[d,ci[d]]
        end
      end
    end
    k = Gridap.Polynomials._set_gradient!(v,s,k,Gridap.TensorValues.SymTensorValue{D})
  end

end



function Gridap.Polynomials._set_gradient!(
  v::AbstractVector{G},s,k,::Type{<:Gridap.TensorValues.SymTensorValue{D}}) where {G,D}
  m = zero(Gridap.TensorValues.Mutable(G))
  # Without the following line modification,
  # it asserts on LoadError: AssertionError: L == (D * (D + 1)) ÷ 2
  # To fix this in Gridap.TensorValues.SymTensorValue
  w = zero(Gridap.TensorValues.SymTensorValue{D,eltype(s),D*(D + 1)÷2})
  z = zero(eltype(s))

  ciw=CartesianIndices(w)
  for c in 1:size(ciw,2) # Go over cols
    for r in 1:c         # Go over upper triangle, current col
      for i in CartesianIndices(m)
        @inbounds m[i] = z
      end
      for i in CartesianIndices(s)
        @inbounds m[i,r,c] = s[i]
        if (r!=c)
          @inbounds m[i,c,r] = s[i]
        end
      end
      @inbounds v[k] = m
      k += 1
    end
  end
  k
end

function Gridap.ReferenceFEs._generate_dof_layout_node_major(
  ::Type{T},
  nnodes::Integer) where T <: Gridap.TensorValues.SymTensorValue
  ncomps = num_components(T)
  V = Gridap.TensorValues.change_eltype(T,Int)
  ndofs = ncomps*nnodes
  dof_to_comp = zeros(Int,ndofs)
  dof_to_node = zeros(Int,ndofs)
  node_and_comp_to_dof = zeros(V,nnodes)
  m = Vector{eltype(V)}(undef,length(V)) # Allocation!
  for node in 1:nnodes
    for comp in 1:ncomps
      o = nnodes*(comp-1)
      dof = node+o
      dof_to_comp[dof] = comp
      dof_to_node[dof] = node
      m[comp] = dof
    end
    node_and_comp_to_dof[node] = Tuple(m)
  end
  (dof_to_node, dof_to_comp, node_and_comp_to_dof)
end

function Gridap.ReferenceFEs._evaluate_lagr_dof!(
  c::AbstractMatrix,
  node_pdof_comp_to_val::AbstractMatrix{<:Gridap.TensorValues.SymTensorValue},
  node_and_comp_to_dof::AbstractVector{<:Gridap.TensorValues.SymTensorValue},
  ndofs,
  ncomps)
  _, npdofs = size(node_pdof_comp_to_val)
  Gridap.Arrays.setsize!(c,(ndofs,npdofs))
  r = c.array
  for node in LinearIndices(node_and_comp_to_dof) # [1,2,3]
    comp_to_dof = node_and_comp_to_dof[node]
    for pdof in 1:npdofs # 9
      comp_to_val = node_pdof_comp_to_val[node,pdof]
      for comp in 1:ncomps # 3
        dof = comp_to_dof.data[comp]
        val = comp_to_val.data[comp]
        r[dof,pdof] = val
      end
    end
  end
  r
end
