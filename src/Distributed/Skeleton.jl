
function SkeletonTriangulation(model::GridapDistributed.DistributedDiscreteModel)
  SkeletonTriangulation(no_ghost,model::GridapDistributed.DistributedDiscreteModel)
end


function _sign_flips(model::GridapDistributed.DistributedDiscreteModel)
  cell_reffes = map_parts(model.models) do m
    basis,reffe_args,reffe_kwargs = ReferenceFE(raviart_thomas,Float64,0)
    reffe = ReferenceFE(m,basis,reffe_args...;reffe_kwargs...)
  end
  GridapDistributed._generate_sign_flips(model,cell_reffes)
end

function SkeletonTriangulation(portion::GridapDistributed.WithGhost,
                               model::GridapDistributed.DistributedDiscreteModel)
   dtrians=map_parts(model.models,_sign_flips(model)) do model, sign_flip
      Skeleton(model,sign_flip)
   end
   GridapDistributed.DistributedTriangulation(dtrians,model)
end

function SkeletonTriangulation(portion::GridapDistributed.NoGhost,
                               model::GridapDistributed.DistributedDiscreteModel)
  cell_gids=get_cell_gids(model)
  dtrians=map_parts(model.models,
                    _sign_flips(model),
                    cell_gids.partition) do model, sign_flip, partition
    cell_to_parent_cell=
      findall([partition.lid_to_part[cell]==partition.part
                 for cell=1:length(partition.lid_to_part)])
    Skeleton(model,cell_to_parent_cell,sign_flip)
  end
  GridapDistributed.DistributedTriangulation(dtrians,model)
end

function get_cell_normal_vector(dskeleton::GridapDistributed.DistributedTriangulation)
  fields=map_parts(dskeleton.trians) do skeleton
    get_cell_normal_vector(skeleton)
  end
  GridapDistributed.DistributedCellField(fields)
end

function get_cell_owner_normal_vector(dskeleton::GridapDistributed.DistributedTriangulation)
  fields=map_parts(dskeleton.trians) do skeleton
    get_cell_owner_normal_vector(skeleton)
  end
  GridapDistributed.DistributedCellField(fields)
end
