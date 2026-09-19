# Dual-axis MEMS micromirror

This notebook builds an invariant-manifold ROM of one vibration mode of a
dual-axis MEMS micromirror and reads its **backbone curve**, the
amplitude-dependent frequency, off the reduced dynamics. It follows the same
public MORFE workflow as [From a mesh to a ROM](../from_a_mesh_to_a_rom/README.md),
on a real device geometry instead of a beam.

The micromirror is a circular plate suspended in three concentric rings. Torsion
bars along x and y join each ring to the next, and the outer ring to the frame,
so the plate can tilt about both axes. The mesh
[`micromirror_device_layer.msh`](micromirror_device_layer.msh) covers the device
layer, 45 thick, with one layer of linear hexahedra; MORFE puts a quadratic
displacement field on it. The device layer is bonded to a handle wafer outside
the backside cavity, radius 1450; the handle itself is not modelled.

Every boundary surface of the mesh is a numbered facet group `surface_<id>`, so
the clamp is chosen in the notebook rather than baked into the mesh.

## Workflow

The notebook stops twice for a look in ParaView:

1. **Clamp.** Section 1 writes the boundary surfaces with their numbers to
   `results/paraview/`. View them with

   ```bash
   paraview --script=show_surfaces.py
   ```

   from this folder, which writes each `surface_id` on its surface (on macOS
   the executable is `/Applications/ParaView-<version>.app/Contents/MacOS/paraview`).
   *Hover Cells On* in the toolbar also shows the number under the cursor. Put
   the numbers in `clamped` and re-run the cell. The default, surface 8, is the
   face bonded to the handle.
2. **Master mode.** Section 3 prints the first ten modes with the tilt of the
   mirror plate about x and y, and writes them to `results/paraview/modes.pvd`.
   Open it, apply *Warp By Vector* with `u`, and step through the modes with the
   time slider: its value is the mode number. Put that number in `p` in
   section 4.

The rest runs as in the beam example: `build_model` with `master = [p]`,
`parametrise`, the backbone from `normal_form_branch`, and `save_rom`. The
observable is the tilt of the mirror plate, the slope of a plane fitted to the
out-of-plane displacement of its face.

```julia
case = SVK.mechanical_model(grid; material = SVK.CubicCrystal(...), ...)
sp = SVK.spectrum(case; nev)
(; model, spectral, meta) = build_model(case; master = [p], spectrum = sp,
    expansion_order = order)
W, R = parametrise(model, spectral, order;
    resonance = ResonanceConfig(style = :complex_normal_form, tol_relative = 0.05))
θ = observable_polynomial(W, θx)   # mirror tilt, a linear functional of W
```

## Run

From the repository root, set up the shared environment once:

```bash
julia setup.jl
```

Then open [`dual_axis_mems_micromirror.ipynb`](dual_axis_mems_micromirror.ipynb)
with the `MORFEExamples` kernel. [`paraview.jl`](paraview.jl) holds the two
ParaView writers the notebook uses.

## Files

| File | Content |
|---|---|
| `micromirror_device_layer.msh` | Hexahedral mesh of the device layer, one facet group per boundary surface |
| `dual_axis_mems_micromirror.ipynb` | Clamp, modes, ROM, backbone |
| `paraview.jl` | Writers for the numbered surfaces and the modes |
| `show_surfaces.py` | ParaView script that labels each surface with its number |
| `results/` | ROM, backbone and ParaView output; `results/paraview/` is not tracked |
