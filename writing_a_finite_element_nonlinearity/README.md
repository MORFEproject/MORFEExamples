# Writing a finite element nonlinearity

Tutorial page: [Writing a finite element nonlinearity](https://morfeproject.github.io/tutorials/fe_nonlinearity.html).

This notebook writes the geometric nonlinearity of the clamped–clamped St.
Venant–Kirchhoff beam of [`from_a_mesh_to_a_rom`](../from_a_mesh_to_a_rom/)
by hand, with Ferrite.jl, in two ways:

1. **As a closure** that assembles a global vector. This is the quickest way to plug any
   finite element code into MORFE, but MORFE calls the closure once per factorisation of
   every monomial, and every call traverses the whole mesh.
2. **Through the element-level interface**: a `FEMMultilinearMap` subtype that describes
   one cell at a time. MORFE then runs the mesh loop itself and batches every evaluation
   of a monomial into one traversal.

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

## Benchmark: batching on larger meshes

`benchmark.jl` measures the batching that the MORFE.jl paper reports, on three generated
meshes, all 1000 × 10 × 24 with quadratic Lagrange displacements:

| Mesh | Elements | Free DOFs |
| ---- | -------- | --------- |
| medium | 60 × 4 × 2 | 16,065 |
| large | 120 × 6 × 3 | 65,247 |
| extra-large | 160 × 8 × 4 | 146,421 |

Both routes run MORFEFerrite's production SVK kernel. The batched route is the built-in
term. The closure route wraps the same term in a closure that traverses the mesh once per
call, with the same factorisations, so the two differ only in the batching.

Nothing is solved. Assembling the right-hand side of a monomial costs the same whatever the
values of the coefficients, so the script fills `W` with random coefficients and times
MORFE's `compute_multilinear_terms!` directly. It times every monomial that a master-mode-1
solve up to order 9 would assemble, in the solver's order. It fails if the two routes
assemble different vectors.

Run it on an otherwise idle machine: no other Julia processes, builds or browsers. Without
the eigensolve and the linear solves it needs about 4 GB of memory.

```bash
cd writing_a_finite_element_nonlinearity
julia benchmark.jl          # writes results/data/batching_benchmark.csv
julia plot_benchmark.jl     # draws it and exports the paper's data
```

`MORFE_BENCH_SMOKE=1 julia benchmark.jl` checks the pipeline in a few minutes, at order 5
with one repeat. Its CSV is marked as a smoke run and is not for publication.

`batching_benchmark.csv` has one row per mesh, route and monomial, in solve order:

- the mesh traversals;
- the right-hand-side time, the fastest of three repeats;
- the bytes allocated during the assembly.

Its first line records the date, the Julia version, the CPU, the number of threads and the
package versions.

`plot_benchmark.jl` plots, against the number of solved monomials, with dashed lines at the
degree boundaries:

- the cumulative time and allocated memory of both routes;
- the cumulative mesh traversals.

It writes `results/figures/batching_benchmark.png`, `results/data/batching_<mesh>.dat` and
`results/data/batching_degrees.tex`. Copy the `.dat` and `.tex` files to the paper's
`figures/data/` folder.
