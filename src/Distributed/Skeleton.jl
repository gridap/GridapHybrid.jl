
function SkeletonTriangulation(model::GridapDistributed.DistributedDiscreteModel)
  SkeletonTriangulation(no_ghost,model::GridapDistributed.DistributedDiscreteModel)
end

function SkeletonTriangulation(portion::GridapDistributed.WithGhost,
                               model::GridapDistributed.DistributedDiscreteModel)
   dtrians=map_parts(model.models) do model
      Skeleton(model)
   end
   GridapDistributed.DistributedTriangulation(dtrians,model)
end

function SkeletonTriangulation(portion::GridapDistributed.NoGhost,
                               model::GridapDistributed.DistributedDiscreteModel)
  dtrians=map_parts(model.models, model.gids.partition) do model, partition
    #cell_to_parent_cell=
    #  findall([partition.lid_to_part[cell]==partition.part
    #             for cell=1:length(partition.lid_to_part)])
    #modelp=DiscreteModelPortion(model, cell_to_parent_cell)
    Skeleton(model)
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
