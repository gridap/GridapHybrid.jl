struct DensifyInnerMostBlockLevel <: Gridap.Fields.Map
end

function Gridap.Arrays.return_cache(k::DensifyInnerMostBlockLevel,
                                    a::Gridap.Fields.ArrayBlock{<:Array{T,M},N}) where {T,M,N}
      @assert all(a.touched)
      block_size=size(a)
      max_M_N=max(M,N)
      densified_size=Vector{Int}(undef,max_M_N)
      current_index=Vector{Int}(undef,N)
      for D=1:N
         current_index[1:D-1].=1
         current_index[D+1:N].=1
         densified_size[D]=0
         for I=1:block_size[D]
           current_index[D]=I
           if (D<=M)
             densified_size[D]+=size(a.array[current_index...])[D]
           else
             densified_size[D]+=1
           end
         end
      end
      if (M>N)
        current_index.=1
        for I=N+1:M
          densified_size[I]=size(a.array[current_index...])[I]
        end
      end
      Gridap.Arrays.CachedArray(Array{T,max_M_N}(undef,Tuple(densified_size)))
end

function Gridap.Arrays.return_value(cache,
                                    k::DensifyInnerMostBlockLevel,
                                    a::Gridap.Fields.ArrayBlock{<:Array{T,M},N}) where {T,M,N}
      @assert all(a.touched)
      max_M_N=max(M,N)
      block_size=size(a)
      scalar_size=size(cache)
      output=cache.array
      current_block_ranges    = Vector{UnitRange{Int}}(undef,max_M_N)
      upper_left_entry_index  = Vector{Int}(undef,max_M_N)
      upper_left_entry_index .= 1
      cinds=CartesianIndices(size(a))
      # if (M>N)
      #   cind=cinds[1]
      #   current_block_size=size(a.array[cind...])
      #   for I=N+1:M
      #     current_block_ranges[I]=1:current_block_size[I]
      #   end
      # end
      for cind in cinds
        current_block_size=size(a.array[cind])
        for (p,s) in enumerate(current_block_size)
            current_block_ranges[p]=upper_left_entry_index[p]:upper_left_entry_index[p]+s-1
            upper_left_entry_index[p]=upper_left_entry_index[p]+s
            if (upper_left_entry_index[p]>scalar_size[p])
              upper_left_entry_index[p]=1
            end
        end
        if (M<N)
          for p=M+1:N
            current_block_ranges[p]=upper_left_entry_index[p]:upper_left_entry_index[p]
            upper_left_entry_index[p]=upper_left_entry_index[p]+1
            if (upper_left_entry_index[p]>scalar_size[p])
              upper_left_entry_index[p]=1
            end
          end
        end
        output[current_block_ranges...]=a.array[cind]
      end
      output
end
