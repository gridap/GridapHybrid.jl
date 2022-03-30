function _merge_bulk_and_skeleton_contributions(
    matcontribs::GridapDistributed.DistributedDomainContribution,
    veccontribs::GridapDistributed.DistributedDomainContribution)
  dmat,dvec=map_parts(matcontribs.contribs,veccontribs.contribs) do matseq, vecseq
    _merge_bulk_and_skeleton_contributions(matseq,vecseq)
  end
  GridapDistributed.DistributedDomainContribution(dmat),
          GridapDistributed.DistributedDomainContribution(dvec)
end

function Gridap.FESpaces._pair_contribution_when_possible(
    obiform::GridapDistributed.DistributedDomainContribution,
    oliform::GridapDistributed.DistributedDomainContribution)
  dmv,dm,dv=map_parts(obiform.contribs,oliform.contribs) do obiform, oliform
    Gridap.FESpaces._pair_contribution_when_possible(obiform,oliform)
  end
  GridapDistributed.DistributedDomainContribution(dmv),
    GridapDistributed.DistributedDomainContribution(dm),
      GridapDistributed.DistributedDomainContribution(dv)
end

function Gridap.FESpaces._collect_cell_matrix_and_vector(
  trial  :: GridapDistributed.DistributedFESpace,
  test   :: GridapDistributed.DistributedFESpace,
  matvec :: GridapDistributed.DistributedDomainContribution,
  mat    :: GridapDistributed.DistributedDomainContribution,
  vec    :: GridapDistributed.DistributedDomainContribution)
  map_parts(GridapDistributed.local_views(trial),
            GridapDistributed.local_views(test),
            matvec.contribs,mat.contribs,vec.contribs) do trial, test, matvec, mat, vec
     Gridap.FESpaces._collect_cell_matrix_and_vector(trial,test,matvec,mat,vec)
  end
end

function _add_static_condensation(
   matvec::GridapDistributed.DistributedDomainContribution,
   bulk_fields, skeleton_fields)
   dmv=map_parts(matvec.contribs) do matvec
    _add_static_condensation(matvec,bulk_fields,skeleton_fields)
  end
  GridapDistributed.DistributedDomainContribution(dmv)
end

function Gridap.FESpaces._attach_dirichlet(
  matvec:: GridapDistributed.DistributedDomainContribution,
  mat   :: GridapDistributed.DistributedDomainContribution,
  uhd)
  dmv,dm=map_parts(matvec.contribs,mat.contribs,uhd.fields) do matvec, mat, uhd
    m,v=Gridap.FESpaces._attach_dirichlet(matvec,mat, uhd)
    m,v
  end
  GridapDistributed.DistributedDomainContribution(dmv),
      GridapDistributed.DistributedDomainContribution(dm)
end


function _compute_hybridizable_from_skeleton_free_dof_values(
    skeleton_fe_function :: GridapDistributed.DistributedCellField,
    trial_hybridizable   :: GridapDistributed.DistributedFESpace,
    test_hybridizable    :: GridapDistributed.DistributedFESpace,
    trial_skeleton       :: GridapDistributed.DistributedFESpace,
    matvec               :: GridapDistributed.DistributedDomainContribution,
    bulk_fields, skeleton_fields)

    function f(skel_fe, trial_hyb, test_hyb, trial_skel, matvec)
      _compute_hybridizable_from_skeleton_free_dof_values(skel_fe,
                                                          trial_hyb,
                                                          test_hyb,
                                                          trial_skel,
                                                          matvec,
                                                          bulk_fields,
                                                          skeleton_fields)
    end
    # f(skeleton_fe_function.fields.parts[1],
    #   trial_hybridizable.part_fe_space.parts[1],
    #   test_hybridizable.part_fe_space.parts[1],
    #   trial_skeleton.spaces.parts[1],
    #   matvec.contribs.parts[1])
    values=map_parts(f,GridapDistributed.local_views(skeleton_fe_function),
             GridapDistributed.local_views(trial_hybridizable),
             GridapDistributed.local_views(test_hybridizable),
             GridapDistributed.local_views(trial_skeleton),
             GridapDistributed.local_views(matvec))
    PVector(values,trial_hybridizable.gids)
end
