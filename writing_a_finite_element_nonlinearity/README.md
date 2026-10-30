# Writing a finite element nonlinearity

Tutorial page: [Writing a finite element nonlinearity](https://morfeproject.github.io/tutorials/fe_nonlinearity.html).

This notebook writes the geometric nonlinearity of the clamped–clamped St.
Venant–Kirchhoff beam of [`from_a_mesh_to_a_rom`](../from_a_mesh_to_a_rom/)
by hand, with Ferrite.jl, in two ways:

1. **As a closure** that assembles a global vector. This is the quickest way to plug any
   finite element code into MORFE, but MORFE calls the closure once per factorisation of
   every monomial, and every call traverses the whole mesh.
2. **Through the element-level interface**: a `FEMMultilinearMap` subtype that describes
   one cell at a time. MORFE then runs the mesh loop itself and assembles every evaluation
   of a monomial in a single traversal. The term opts into flux summation (`fem_flux_zero`,
   `flux_qp`, `contract_qp!`), so at each quadrature point MORFE sums the stresses of all
   factorisations of both terms and tests them once against the shape functions.

The physics is the same few lines in both routes. They are the quadratic and cubic parts
of the first Piola–Kirchhoff stress, `P2(a, b)` and `P3(a, b, c)`, polarised so that they
accept distinct gradients:

```julia
σ(E) = λ * tr(E) * one(E) + 2μ * E
ε(a) = symmetric(a)
E2(a, b) = symmetric(transpose(a) ⋅ b + transpose(b) ⋅ a) / 4
P2(a, b) = -(σ(E2(a, b)) + (a ⋅ σ(ε(b)) + b ⋅ σ(ε(a))) / 2)
P3(a, b, c) = -(a ⋅ σ(E2(b, c)) + b ⋅ σ(E2(a, c)) + c ⋅ σ(E2(a, b))) / 3
```

The mesh, the linear matrices and the spectral data come from `SVK.mechanical_model` and
`build_model`, so the notebook can check both routes against the backend's built-in term.
All three ROMs agree to round-off. At order 9 the closure traverses the mesh 486 times and
the element-level term 28 times, once per solved monomial.

## Run

From the repository root, set up the shared environment once:

```bash
julia setup.jl
```

Open and execute
[`writing_a_finite_element_nonlinearity.ipynb`](writing_a_finite_element_nonlinearity.ipynb)
with the `MORFEExamples` kernel selected. It reads the mesh from
`../from_a_mesh_to_a_rom/clamped_clamped_beam.msh` and runs in a few minutes.

From a shell:

```bash
cd writing_a_finite_element_nonlinearity
python3 -m nbconvert --execute --to notebook --inplace writing_a_finite_element_nonlinearity.ipynb
```

The notebook's last cell writes `results/data/rhs_cost.csv`: the mesh traversals and the
right-hand-side time of both routes for orders 3 to 9, the fastest of three runs each. The
tutorial page draws it.

## Benchmark: single-pass assembly on larger meshes

`benchmark.jl` measures the single-pass assembly that the MORFE.jl paper reports. By default it uses the
medium mesh, which is what the paper shows; the larger ones are optional
(`MORFE_BENCH_MESHES=medium,large,extra-large`). The small one has the 4977 DOFs of the mesh of
`from_a_mesh_to_a_rom`, and gives the tutorial page the costs of MORFEFerrite's kernel:

```bash
MORFE_BENCH_MESHES=small MORFE_BENCH_ROMS=2 MORFE_BENCH_ORDERS="2=>9" \
    MORFE_BENCH_OUT=benchmark.v2.csv julia benchmark.jl
```

All four are generated, 1000 × 10 × 24, with quadratic Lagrange displacements:

| Mesh | Elements | Free DOFs |
| ---- | -------- | --------- |
| small | 40 × 3 × 1 | 4,977 |
| medium | 60 × 4 × 2 | 16,065 |
| large | 120 × 6 × 3 | 65,247 |
| extra-large | 160 × 8 × 4 | 146,421 |

It also uses four reduced dimensions. Each goes to the order whose multiindex set holds about
200 monomials:

| n | Reduced coordinates | Order | Monomials | Assembled | Closure traversals |
| - | ------------------- | ----- | --------- | --------- | ------------------ |
| 2 | one mode pair | 19 | 209 | 108 | 16,003 |
| 3 | one mode pair and a parameter | 9 | 219 | 122 | 4,778 |
| 4 | two mode pairs | 6 | 209 | 107 | 1,697 |
| 5 | two mode pairs and a parameter | 5 | 251 | 132 | 1,409 |

Both routes run MORFEFerrite's production SVK kernel. The single-pass route is the built-in
term. The closure route wraps the same term in a closure that traverses the mesh once per
call, with the same factorisations, so the two differ only in the single pass.

Nothing is solved. Assembling the right-hand side of a monomial costs the same whatever the
values of the coefficients, so the script works as follows:

- For each mesh and n it allocates one `W`, fills it with random vectors for all its
  monomials at once, and refills it in place for each of three random realisations.
- It times MORFE's `compute_multilinear_terms!` for every monomial the solver would assemble,
  in the solver's order, and averages time and allocations over the realisations.
- Each (mesh, n) is a single run to its order. The solve is graded, so that run contains
  every lower order.
- It fails if the two routes assemble different vectors.

Run it on an otherwise idle machine. On an Apple M2 Pro the medium mesh takes about 40
minutes; one closure traversal takes 0.015 s there, against 0.075 s on the large mesh and
0.16 s on the extra-large one. All three meshes take about 8 hours, and the extra-large one
needs about 2 GB of memory. The ratios the plots show change by only a few per cent from the
medium to the larger meshes, so the medium mesh suffices.

```bash
cd writing_a_finite_element_nonlinearity
julia benchmark.jl          # writes results/data/single_pass_benchmark.csv
julia plot_benchmark.jl     # draws it and exports the paper's data
```

Overrides:

- `MORFE_BENCH_MESHES=medium`
- `MORFE_BENCH_ROMS=2,3`
- `MORFE_BENCH_OUT=<csv path>`, to write elsewhere than `results/data/single_pass_benchmark.csv`
- `MORFE_BENCH_ORDERS="2=>19,3=>9,4=>6,5=>5"`
- `MORFE_BENCH_REALISATIONS=3`
- `MORFE_BENCH_PROJECT=<environment>`, to benchmark local checkouts added with `Pkg.develop`

`single_pass_benchmark.csv` has one row per mesh, n, route (`closure` or `single_pass`) and
monomial, in solve order:

- the mesh traversals;
- the mean assembly time;
- the mean bytes allocated during the assembly;
- the workspace the route keeps for the assembly, which is its factorisation cache and
  buffers (not the FEM terms and the mesh that the single-pass cache refers to).

Its first line records the date, the Julia version, the CPU, the number of threads and the
package versions.

`plot_benchmark.jl` plots four panels, in a 2 × 2 grid, against the number of monomials of the expansion,
in the solver's order. The solver assembles a right-hand side for one monomial of each
complex-conjugate pair and none for the linear monomials, so the cumulative curves stay
flat across those. The panels are:

- the speed-up of the single-pass route, closure time / single-pass time, accumulated;
- the bytes the closure allocates during assembly, accumulated, on a log scale (the single
  pass allocates none; its total, if any, is in the axis label);
- the order reached, the same for every mesh, on a log scale;
- the factorisations evaluated (one closure call each), the same for every mesh, on a log
  scale with a gridline at each decade.

The monomial axis is logarithmic, from 1, labelled at the decades with tick marks at every
1–9 multiple. The log axes of the memory
and factorisation panels are labelled at the decades too. Colour is n, and the medium mesh is drawn solid; a run with more meshes adds
dashed and dotted curves for the large and extra-large ones. Markers end each order and cycle through circle, star, triangle, cross and
square. In the speed-up, memory and factorisation panels, faint grey lines join the same
order across n.

It writes these files:

- `results/figures/single_pass_benchmark.png`;
- one `results/data/single_pass_assembly_<mesh>_n<n>.dat` per mesh and n;
- one `results/data/single_pass_assembly_orders_n<n>.dat` per n, holding the degree boundaries;
- one `results/data/single_pass_assembly_<mesh>_isoorders.dat` per mesh, holding the
  same-order markers across n.

Copy the `.dat` files to the paper's `figures/data/` folder; the names already match.

`phase_breakdown.jl` splits the time of both routes, over the same expansions, into five
phases per element: setup (Jacobian and shape-function gradients), gradients, stresses,
contraction against the test functions and assembly. Timers inside the loops would distort
phases that last a fraction of a microsecond, so it runs each route five times, stopping
after one more phase each time, and takes the differences; the full run is checked against
MORFE's own right-hand side. It takes about 20 minutes and writes
`results/data/single_pass_phases.csv`, from which the paper's phase bars are drawn.
