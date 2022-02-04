
# These codes should go eventually to Gridap at the end!!!


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
