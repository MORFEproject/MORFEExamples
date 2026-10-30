# Right-hand-side assembly cost of the two finite element routes on the clamped-clamped
# beam: a closure that traverses the mesh once per evaluation, that is once per
# factorisation, against MORFE's single-pass assembly through the element-level interface,
# which traverses it once per monomial. Both routes run MORFEFerrite's production SVK
# kernel, so the comparison isolates the single pass. Three mesh sizes and four reduced dimensions n:
#   n = 2  one mode pair              n = 4  two mode pairs
#   n = 3  one mode pair + parameter  n = 5  two mode pairs + parameter
#
# Nothing is solved. Assembling the right-hand side of a monomial costs the same whatever the
# values of the coefficients, so for every mesh and n one W is allocated, filled with random
# vectors for all its monomials at once, and refilled in place for each of the random
# realisations. MORFE's `compute_multilinear_terms!` is timed for every monomial the solver
# would assemble, in the solver's order; time and allocations are averaged over the
# realisations.
#
# Each (mesh, n) is a single run to the order whose multiindex set holds about 200 monomials
# (ORDERS below). The solve is graded, so that run contains every lower order: its
# cumulative curves read off any of them at the degree boundaries.
#
#   julia benchmark.jl    # on an otherwise idle machine
#
# Optional overrides: MORFE_BENCH_ORDERS="2=>19,3=>9,4=>6,5=>5"  MORFE_BENCH_REALISATIONS=3
# MORFE_BENCH_MESHES=small,medium,large,extra-large (default: medium)  MORFE_BENCH_ROMS=2,3,
# MORFE_BENCH_OUT=<csv path>, and MORFE_BENCH_PROJECT=<environment> to benchmark, for instance,
# local MORFE and MORFEFerrite checkouts added with Pkg.develop. By default the script activates the repository's shared
# environment itself; run `julia setup.jl` from the repository root once beforehand.
#
# Writes results/data/single_pass_benchmark.csv, one row per (mesh, n, route, monomial):
# the mesh traversals, the mean assembly time and allocated bytes, and the workspace the
# route keeps for the assembly. plot_benchmark.jl draws it.
using Pkg: Pkg
Pkg.activate(get(ENV, "MORFE_BENCH_PROJECT", dirname(@__DIR__)))
using MORFE, MORFEFerrite, Ferrite, Dates, Printf, Random
using MORFE.Polynomials: DensePolynomial
using MORFE.Multiindices: find_in_set
using MORFE.MultilinearTerms: build_multilinear_terms_cache, compute_multilinear_terms!
const SVK = MORFEFerrite.StructuralSVK

# Elements along the length (1000), the thickness (10) and the width (24) of the beam, with
# quadratic Lagrange displacements; small, 40 x 3 x 1, has the 4977 dofs of the mesh of
# from_a_mesh_to_a_rom, on which the tutorial page draws this kernel's costs.
const ALL_MESHES = ["small" => (40, 3, 1), "medium" => (60, 4, 2), "large" => (120, 6, 3),
    "extra-large" => (160, 8, 4)]
# Reduced dimension n => expansion order, chosen so that every multiindex set holds about
# 200 monomials: 209, 219, 209 and 251 (of which the solver assembles 108, 122, 107 and 132;
# the closure then traverses the mesh 16003, 4778, 1697 and 1409 times).
const ORDERS = Dict(parse(Int, a) => parse(Int, b)
for (a, b) in split.(split(get(ENV, "MORFE_BENCH_ORDERS", "2=>19,3=>9,4=>6,5=>5"), ","), "=>"))
const REALISATIONS = parse(Int, get(ENV, "MORFE_BENCH_REALISATIONS", "3"))
const ROMS = parse.(Int, split(get(ENV, "MORFE_BENCH_ROMS", "2,3,4,5"), ","))
# The paper uses the medium mesh: the ratios it reports hardly change on the larger ones.
const MESHES = let names = split(get(ENV, "MORFE_BENCH_MESHES", "medium"), ",")
    filter(m -> first(m) in names, ALL_MESHES)
end

material = SVK.SVKMaterial(E = 160e3, ν = 0.22, ρ = 2.32e-3)
damping = SVK.RayleighDamping(α = 0.0, β = 0.0)

function beam(n)
    grid = generate_grid(Hexahedron, n, Vec((0.0, 0.0, 0.0)), Vec((1000.0, 10.0, 24.0)))
    addfacetset!(grid, "Dirichlet", x -> abs(x[1]) < 1e-8 || abs(x[1] - 1000) < 1e-8)
    return SVK.mechanical_model(grid; material, damping, dirichlet = "Dirichlet")
end

# The production term behind a closure: every call is one full mesh traversal. The closure
# is tuned as far as one evaluation allows, so the comparison does not flatter the single pass:
#   - a repeated argument, as in T(a, a, b), reuses the gradients of its first occurrence;
#   - per element, the stresses of all quadrature points are formed first (`flux_qp`) and
#     tested against the shape functions afterwards (`contract_qp!`). Point by point, as
#     `accumulate_qp!` does it, each test waits for its stress and the kernel runs 1.4x
#     slower; the arithmetic is the same.
const passes = Ref(0)

struct PerFactorisation{T, B, P, F}
    t::T
    ∇u::B
    P::P
    Fe::F
end
PerFactorisation(t) = PerFactorisation(t, similar(fem_qp_buffer(t), t.deg, fem_n_qp(t)),
    Vector{typeof(fem_flux_zero(t))}(undef, fem_n_qp(t)),
    zeros(ComplexF64, fem_ndofs_per_cell(t)))

function (u::PerFactorisation)(res, args::Vararg{Any, N}) where {N}
    passes[] += 1
    t = u.t
    n_qp = fem_n_qp(t)
    for element in fem_elements(t)
        fem_reinit!(element, t)
        for k in 1:N
            j = 0
            for i in 1:(k - 1)
                args[i] === args[k] && (j = i; break)
            end
            if j == 0
                scatter_qp!(view(u.∇u, k, :), args[k], element, t)
            else
                copyto!(view(u.∇u, k, :), view(u.∇u, j, :))
            end
        end
        fill!(u.Fe, 0)
        for q in 1:n_qp
            u.P[q] = flux_qp(ntuple(k -> u.∇u[k, q], Val(N)), element, q, t)
        end
        for q in 1:n_qp
            contract_qp!(u.Fe, u.P[q], element, q, fem_getdetJdV(element, q, t), t)
        end
        assemble_element!(res, u.Fe, element, t)
    end
    return res
end

closure(terms) = map(t -> MultilinearMap(PerFactorisation(t); multiindex = t.multiindex,
        multiplicity_external = 0, fully_asymmetric = false), terms)

# Complex-conjugate partner of an exponent: the indices of each mode pair swap, a parameter
# maps to itself. Variables: (z₁, z̄₁[, z₂, z̄₂][, η]).
partner(e, n) = Tuple(n >= 4 ? (e[2], e[1], e[4], e[3], e[5:end]...) : (e[2], e[1], e[3:end]...))

# The monomials the solver assembles: degree 2 and up, and of each conjugate pair only the
# first in the solve order, since the solver fills its partner by conjugation.
assembled(mset, n) = [l for l in eachindex(mset.exponents)
                      if sum(mset.exponents[l]) >= 2 &&
                         find_in_set(mset, partner(Tuple(mset.exponents[l]), n)) >= l]

function assemble!(buf, m, l, W, cache)
    before = passes[]
    r = @timed compute_multilinear_terms!(buf, m, l, W, cache)
    return (r.time, r.bytes, passes[] - before)
end

# Workspace a route keeps for the assembly: its factorisation cache and buffers, plus, for
# the closure, the per-term gradient, stress and element buffers, and, for the single pass,
# the per-coefficient strains of the SVK terms' joint stress. The single-pass cache also
# refers to the FEM terms themselves (the `terms` of each group in `fem_groups`); they, and
# the mesh they hold, are left out. (`summarysize`'s `exclude` cannot do it: the terms sit
# inline in those tuples.)
function workspace(cache, terms)
    closure_buffers(t) = Base.summarysize(t.f!.∇u) + Base.summarysize(t.f!.P) +
                         Base.summarysize(t.f!.Fe)
    term_buffers(t) = t isa MultilinearMap && t.f! isa PerFactorisation ? closure_buffers(t) :
                      t isa FEMMultilinearMap && hasproperty(t, :strains) ?
                      Base.summarysize(t.strains) : 0
    fields_size(x, skip) = sum(f -> f === skip ? 0 : Base.summarysize(getfield(x, f)),
        fieldnames(typeof(x)))
    group_size(groups) = sum(g -> fields_size(g, :terms), groups; init = 0)
    cache_size = fields_size(cache, :fem_groups) + group_size(cache.fem_groups)
    return cache_size + sum(term_buffers, terms)
end

println("realisations $REALISATIONS, meshes ", join(first.(MESHES), ", "), ", ",
    join(("n = $n to order $(ORDERS[n])" for n in ROMS), ", "))
out = get(ENV, "MORFE_BENCH_OUT",
    joinpath(@__DIR__, "results", "data", "single_pass_benchmark.csv"))
mkpath(dirname(out))
compiled = Set{Int}()
open(out, "w") do io
    println(io, "# $(now()), Julia $(VERSION), $(Sys.cpu_info()[1].model), ",
        "$(Threads.nthreads()) thread(s), MORFE $(pkgversion(MORFE)), ",
        "MORFEFerrite $(pkgversion(MORFEFerrite)), mean of $REALISATIONS realisation(s)")
    println(io, "mesh,elements,n_dofs,n,order,route,k,monomial_idx,degree,exponents,passes,",
        "rhs_time_s,rhs_alloc_bytes,workspace_bytes")
    for (name, elements) in MESHES
        case = beam(elements)
        N = size(case.B[1], 1)
        for n in ROMS
            order = ORDERS[n]
            mset = build_multiindex_set(order, n)
            ls = assembled(mset, n)
            terms = Tuple(case.term_factory(d) for d in (2, 3))
            models = (NthOrderModel(case.B, closure(terms)), NthOrderModel(case.B, terms))
            W = Parametrisation(DensePolynomial(zeros(ComplexF64, N, 2, length(mset)), mset), 0)
            caches = map(m -> build_multilinear_terms_cache(m, W), models)
            bufs = (zeros(ComplexF64, N), zeros(ComplexF64, N))
            if !(n in compiled)   # compile every specialisation: one untimed sweep
                for l in ls, j in 1:2
                    assemble!(bufs[j], models[j], l, W, caches[j])
                end
                push!(compiled, n)
            end
            time = zeros(length(ls), 2)
            bytes = zeros(length(ls), 2)
            calls = zeros(Int, length(ls))
            for r in 1:REALISATIONS
                randn!(Xoshiro(r), W.poly.coefficients)   # refill the same W in place
                GC.gc()
                for (k, l) in enumerate(ls)
                    for j in 1:2
                        t, b, p = assemble!(bufs[j], models[j], l, W, caches[j])
                        time[k, j] += t / REALISATIONS
                        bytes[k, j] += b / REALISATIONS
                        j == 1 && (calls[k] = p)
                    end
                    dev = maximum(abs, bufs[1] - bufs[2]) / maximum(abs, bufs[2])
                    dev < 1e-10 ||
                        error("$name, n = $n, monomial $l: closure and single pass differ, $dev")
                end
            end
            # after the assembly, which grows the single pass's per-coefficient strains
            spaces = map((c, m) -> workspace(c, m.nonlinear_terms), caches, models)
            for (j, route) in enumerate(("closure", "single_pass")), (k, l) in enumerate(ls)
                e = mset.exponents[l]
                println(io, join((name, join(elements, "x"), N, n, order, route, k, l, sum(e),
                        join(e, "_"), j == 1 ? calls[k] : 1, time[k, j], bytes[k, j],
                        spaces[j]), ","))
            end
            flush(io)
            @printf("%-12s %7d dofs  n = %d  order %d  %4d monomials  traversals %6d  rhs %9.3f s / %9.3f s  (x%.2f)  workspace %.1f / %.1f MiB\n",
                name, N, n, order, length(ls), sum(calls), sum(time[:, 1]), sum(time[:, 2]),
                sum(time[:, 1]) / sum(time[:, 2]), spaces[1] / 2^20, spaces[2] / 2^20)
            W = caches = models = nothing
            GC.gc()
        end
    end
end
println("wrote ", out)
