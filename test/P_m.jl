# uh : Trial Basis of L2-projection source space (bulk) or FEFunction
# uhΓ: Trial Basis of L2-projection target space (skeleton)
# μ  : Test  Basis of L2-projection target space (skeleton)
function Pₘ(uh, uhΓ, μ, d∂K)
  A = ∫(μ⋅uhΓ)d∂K
  B = ∫(μ⋅uh)d∂K
  # [c][f][bμ,buhΓ==bμ][f,f]
  A_array = _remove_sum_facets(_remove_densify(Gridap.CellData.get_array(A)))
  # [c][f][bμ,buh][f,1] (if uh Trial Basis)
  # [c][f][bμ]    [f]   (if uh FE Function)
  B_array = _remove_sum_facets(_remove_densify(Gridap.CellData.get_array(B)))
  # [c][f][1,buh][1]   (if uh Trial Basis)
  # [c][f]             (if uh FE Function)
  pm_uh_dofs = lazy_map(GridapHybrid.compute_bulk_to_skeleton_l2_projection_dofs, A_array, B_array)
  # [c][f][1,buhΓ][1,f]
  uhΓ_d∂K = Gridap.CellData.change_domain(uhΓ, d∂K.quad.trian, ReferenceDomain())
  uhΓ_d∂K_data = Gridap.CellData.get_data(uhΓ_d∂K)
  # [c][f][1,buh][1] (if uh Trial Basis)
  # [c][f]           (if uh FE Function)
  field_array = lazy_map(GridapHybrid.setup_bulk_to_skeleton_l2_projected_fields,
    pm_uh_dofs, uhΓ_d∂K_data)
  Gridap.CellData.GenericCellField(field_array, d∂K.quad.trian, ReferenceDomain())
end

function _remove_densify(
     a::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.DensifyInnerMostBlockLevelMap}})
  a.args[1]
end

function _remove_sum_facets(
  a::Gridap.Arrays.LazyArray{<:Fill{<:SumFacetsMap}})
  a.args[1]
end
