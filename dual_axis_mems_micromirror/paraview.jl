# ParaView output for the micromirror notebook: the numbered boundary surfaces
# and the vibration modes. MORFEFerrite's own ParaView writers accept wedges and
# 27-node hexahedra only, so these use Ferrite and WriteVTK directly, which
# handle the linear hexahedra of this mesh.

using WriteVTK: MeshCell, VTKCellData, VTKCellTypes, VTKPointData, paraview_collection,
                vtk_grid

"""
    surface_ids(grid) -> Vector{Int}

Numbers of the facet sets `surface_<id>` of `grid`, sorted.
"""
surface_ids(grid) = sort!([parse(Int, chopprefix(name, "surface_"))
                           for name in keys(Ferrite.getfacetsets(grid))
                           if startswith(name, "surface_")])

"""
    write_surfaces(dir, grid, clamped)

Write the boundary of `grid` to `dir`:

- `surfaces.vtu`: every face of the facet sets `surface_<id>`, with cell data
  `surface_id`, `clamped` (1 on the surfaces listed in `clamped`) and
  `color_key`, a shuffle of the ids that gives neighbouring surfaces
  well-separated colours;
- `surface_labels.vtu`: one point per surface, on the surface, carrying its
  `surface_id`.
"""
function write_surfaces(dir, grid, clamped)
    ids = surface_ids(grid)
    points = reduce(hcat, [collect(Ferrite.get_node_coordinate(n)) for n in Ferrite.getnodes(grid)])
    cells, cell_id, label = MeshCell[], Int[], Int[]
    for s in ids
        nodes = Int[]
        for facet in Ferrite.getfacetset(grid, "surface_$s")
            cell, local_facet = facet[1], facet[2]
            face = collect(Ferrite.facets(Ferrite.getcells(grid, cell))[local_facet])
            push!(cells, MeshCell(VTKCellTypes.VTK_QUAD, face))
            push!(cell_id, s)
            append!(nodes, face)
        end
        centre = sum(points[:, nodes]; dims = 2) ./ length(nodes)
        push!(label, nodes[argmin(vec(sum(abs2, points[:, nodes] .- centre; dims = 1)))])
    end
    vtk_grid(joinpath(dir, "surfaces"), points, cells) do vtk
        vtk["surface_id", VTKCellData()] = cell_id
        vtk["clamped", VTKCellData()] = Int.(cell_id .∈ Ref(clamped))
        vtk["color_key", VTKCellData()] = mod.(cell_id .* sqrt(2), 1.0)
    end
    vtk_grid(joinpath(dir, "surface_labels"), points[:, label],
        [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i in eachindex(ids)]) do vtk
        vtk["surface_id", VTKPointData()] = ids
    end
    return ids
end

"""
    mode_shape(sp, p) -> Vector{Float64}

Displacement of mode pair `p` of the spectrum `sp` on the free degrees of
freedom, made real and scaled so that its largest entry is 1.
"""
function mode_shape(sp, p)
    ϕ = sp.eigenmodes[:, 1, 2p - 1]
    return real.(ϕ ./ ϕ[argmax(abs.(ϕ))])
end

"""
    write_modes(dir, case, sp)

Write every mode pair of `sp` to `dir` as `mode_<p>.vtu`, with point data `u`,
and collect them in `modes.pvd`, whose time value is the mode number `p`.
Constrained degrees of freedom are zero.
"""
function write_modes(dir, case, sp)
    dh = case.info.dh
    pvd = paraview_collection(joinpath(dir, "modes"))
    for p in 1:(length(sp.eigenvalues) ÷ 2)
        ϕ = mode_shape(sp, p)
        u = zeros(Ferrite.ndofs(dh))
        for (dof, row) in case.info.free_to_local
            u[dof] = ϕ[row]
        end
        Ferrite.VTKGridFile(joinpath(dir, "mode_$(lpad(p, 2, '0'))"), dh) do vtk
            Ferrite.write_solution(vtk, dh, u)
            pvd[p] = vtk
        end
    end
    close(pvd)
    return joinpath(dir, "modes.pvd")
end
