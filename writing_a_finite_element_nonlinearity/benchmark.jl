# Right-hand-side assembly cost of the two finite element routes on the clamped-clamped
# beam, at three mesh sizes: a closure that traverses the mesh once per evaluation, that is
# once per factorisation, against MORFE's batched element-level interface, which traverses
# it once per monomial. Both routes run MORFEFerrite's production SVK kernel, so the
# comparison isolates the batching.
#
# Nothing is solved. The cost of assembling the right-hand side of a monomial depends on the
# mesh and on its factorisations, not on the values of the coefficients, so W is filled
# with random coefficients and MORFE's `compute_multilinear_terms!` is timed directly, for
# every monomial the solver would assemble, in the solver's order.
#
#   julia benchmark.jl                      # the full run, on an otherwise idle machine
#   MORFE_BENCH_SMOKE=1 julia benchmark.jl  # pipeline check: order 5, one repeat
#
# Optional overrides: MORFE_BENCH_MESHES=medium,large  MORFE_BENCH_ORDER=9
# MORFE_BENCH_REPEATS=3, and MORFE_BENCH_PROJECT=<environment> to benchmark, for instance,
# local MORFE and MORFEFerrite checkouts added with Pkg.develop. By default the script
# activates the repository's shared environment itself; run `julia setup.jl` from the
# repository root once beforehand.
#
# Writes results/data/batching_benchmark.csv, one row per (mesh, route, monomial): the mesh
# traversals, the assembly time (fastest of the repeats) and the bytes allocated during the
# assembly. plot_benchmark.jl draws it.
using Pkg: Pkg
Pkg.activate(get(ENV, "MORFE_BENCH_PROJECT", dirname(@__DIR__)))
using MORFE, MORFEFerrite, Ferrite, Dates, Printf, Random
using MORFE.Polynomials: DensePolynomial
using MORFE.Multiindices: find_in_set
using MORFE.MultilinearTerms: build_multilinear_terms_cache, compute_multilinear_terms!
const SVK = MORFEFerrite.StructuralSVK

# Elements along the length (1000), the thickness (10) and the width (24) of the beam, with
# quadratic Lagrange displacements; 40 x 3 x 1 is the 4977-dof mesh of from_a_mesh_to_a_rom.
const ALL_MESHES = ["medium" => (60, 4, 2), "large" => (120, 6, 3), "extra-large" => (160, 8, 4)]

const SMOKE = get(ENV, "MORFE_BENCH_SMOKE", "0") == "1"
const ORDER = parse(Int, get(ENV, "MORFE_BENCH_ORDER", SMOKE ? "5" : "9"))
const REPEATS = parse(Int, get(ENV, "MORFE_BENCH_REPEATS", SMOKE ? "1" : "3"))
const MESHES = let names = split(get(ENV, "MORFE_BENCH_MESHES", join(first.(ALL_MESHES), ",")), ",")
    filter(m -> first(m) in names, ALL_MESHES)
end

material = SVK.SVKMaterial(E = 160e3, ν = 0.22, ρ = 2.32e-3)
damping = SVK.RayleighDamping(α = 0.0, β = 0.0)

function beam(n)
    grid = generate_grid(Hexahedron, n, Vec((0.0, 0.0, 0.0)), Vec((1000.0, 10.0, 24.0)))
    addfacetset!(grid, "Dirichlet", x -> abs(x[1]) < 1e-8 || abs(x[1] - 1000) < 1e-8)
    return SVK.mechanical_model(grid; material, damping, dirichlet = "Dirichlet")
end

# The production term behind a closure: every call is one full mesh traversal.
const passes = Ref(0)

struct Unbatched{T, B, F}
    t::T
    ∇u::B
    Fe::F
end
Unbatched(t) = Unbatched(t, similar(fem_qp_buffer(t), t.deg, fem_n_qp(t)),
    zeros(ComplexF64, fem_ndofs_per_cell(t)))

function (u::Unbatched)(res, args::Vararg{Any, N}) where {N}
    passes[] += 1
    t = u.t
    for element in fem_elements(t)
        fem_reinit!(element, t)
        for k in 1:N
            scatter_qp!(view(u.∇u, k, :), args[k], element, t)
        end
        fill!(u.Fe, 0)
        for q in 1:fem_n_qp(t)
            accumulate_qp!(u.Fe, ntuple(k -> u.∇u[k, q], Val(N)), 1, element, q,
                fem_getdetJdV(element, q, t), t)
        end
        assemble_element!(res, u.Fe, element, t)
    end
    return res
end

closure(terms) = map(t -> MultilinearMap(Unbatched(t); multiindex = t.multiindex,
        multiplicity_external = 0, fully_asymmetric = false), terms)

# The monomials the solver assembles: degree 2 and up, and of each complex-conjugate pair
# (z^a zbar^b, z^b zbar^a) of the master mode only the first in the solve order, since the
# solver fills its partner by conjugation.
function assembled(mset)
    return [l for l in eachindex(mset.exponents)
            if sum(mset.exponents[l]) >= 2 &&
               find_in_set(mset, reverse(Tuple(mset.exponents[l]))) >= l]
end

# One assembly, timed: (time, bytes, traversals).
function assemble!(buf, m, l, W, cache)
    before = passes[]
    r = @timed compute_multilinear_terms!(buf, m, l, W, cache)
    return (r.time, r.bytes, passes[] - before)
end

println("orders 2..$ORDER, master mode 1, $REPEATS repeat(s), meshes ",
    join(first.(MESHES), ", "))
out = joinpath(@__DIR__, "results", "data", "batching_benchmark.csv")
mkpath(dirname(out))
open(out, "w") do io
    println(io, "# $(now()), Julia $(VERSION), $(Sys.cpu_info()[1].model), ",
        "$(Threads.nthreads()) thread(s), MORFE $(pkgversion(MORFE)), ",
        "MORFEFerrite $(pkgversion(MORFEFerrite)), order $ORDER, master mode 1, ",
        "fastest of $REPEATS run(s)", SMOKE ? ", SMOKE RUN: not for publication" : "")
    println(io, "mesh,elements,n_dofs,route,k,monomial_idx,degree,exponents,passes,",
        "rhs_time_s,rhs_alloc_bytes")
    for (i, (name, n)) in enumerate(MESHES)
        case = beam(n)
        N = size(case.B[1], 1)
        mset = build_multiindex_set(ORDER, 2)
        terms = Tuple(case.term_factory(d, length(mset)) for d in (2, 3))
        models = ("closure" => NthOrderModel(case.B, closure(terms)),
            "batched" => NthOrderModel(case.B, terms))
        W = Parametrisation(DensePolynomial(randn(Xoshiro(1), ComplexF64, N, 2, length(mset)),
                mset), 0)
        caches = map(((_, m),) -> build_multilinear_terms_cache(m, W), models)
        bufs = (zeros(ComplexF64, N), zeros(ComplexF64, N))
        ls = assembled(mset)
        if i == 1   # compile every specialisation once: one untimed sweep of both routes
            for l in ls, (j, (_, m)) in enumerate(models)
                assemble!(bufs[j], m, l, W, caches[j])
            end
        end
        GC.gc()
        best = [fill((Inf, 0, 0), length(ls)) for _ in models]
        for _ in 1:REPEATS, (k, l) in enumerate(ls), j in eachindex(models)
            r = assemble!(bufs[j], last(models[j]), l, W, caches[j])
            r[1] < best[j][k][1] && (best[j][k] = r)
            if j == 2   # both routes must assemble the same vector
                dev = maximum(abs, bufs[1] - bufs[2]) / maximum(abs, bufs[2])
                dev < 1e-10 || error("$name, monomial $l: closure and batched differ, $dev")
            end
        end
        for (j, (route, _)) in enumerate(models), (k, l) in enumerate(ls)
            t, b, p = best[j][k]
            e = mset.exponents[l]
            println(io, join((name, join(n, "x"), N, route, k, l, sum(e), join(e, "_"),
                    route == "batched" ? 1 : p, t, b), ","))
        end
        flush(io)
        tc, tb = sum(first.(best[1])), sum(first.(best[2]))
        @printf("%-12s %7d dofs  traversals %5d / %3d  rhs %8.3f s / %8.3f s  (x%.2f)\n",
            name, N, sum(last.(best[1])), length(ls), tc, tb, tc / tb)
    end
end
println("wrote ", out)
