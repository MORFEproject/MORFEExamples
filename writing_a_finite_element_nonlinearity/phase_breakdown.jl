# Where the right-hand-side time of each route goes, phase by phase, over the whole expansion
# of benchmark.jl (medium beam mesh, n = 2, 3, 4, 5 at orders 19, 9, 6, 5): the finite element
# closure (M3), which traverses the mesh once per factorisation, against MORFE's single-pass
# assembly (M4), which traverses it once per monomial. The phases, per element:
#
#   setup        the element loop and the cell's reinitialisation (Jacobian, shape-function
#                gradients and integration weights at the quadrature points);
#   gradients    the displacement gradient of each coefficient at every quadrature point;
#   stresses     the polarised first Piola-Kirchhoff stresses at every quadrature point;
#   contraction  the stress tested against the gradients of all shape functions;
#   assembly     the element vector added to the global right-hand side.
#
# Timers inside the loops would distort phases that last a fraction of a microsecond, so the
# phases are measured by cumulative ablation: each route runs at five levels, level k doing
# the first k + 1 phases, and a phase's time is the difference of consecutive levels. Level 4
# is the full route, checked against MORFE's own result. Skipped phases leave their inputs
# in place and every computed quantity is stored or returned, so the compiler drops nothing.
# M3 is benchmark.jl's tuned closure (repeated arguments reuse their gradients; the stresses
# of all quadrature points of an element come first, then their contractions); M4 is MORFE's
# `_replay_fem!` with MORFEFerrite's joint SVK stress, copied phase by phase.
#
#   julia phase_breakdown.jl    # on an otherwise idle machine, about 20 minutes
#
# Overrides as in benchmark.jl: MORFE_BENCH_ORDERS, MORFE_BENCH_ROMS, MORFE_BENCH_PROJECT.
# Writes results/data/single_pass_phases.csv: one row per n, route and level, with the time
# summed over every monomial the solver assembles.
using Pkg: Pkg
Pkg.activate(get(ENV, "MORFE_BENCH_PROJECT", dirname(@__DIR__)))
using MORFE, MORFEFerrite, Ferrite, Dates, Printf, Random
using MORFE.Polynomials: DensePolynomial
using MORFE.Multiindices: find_in_set
using MORFE.MultilinearTerms: build_multilinear_terms_cache, compute_multilinear_terms!
const SVK = MORFEFerrite.StructuralSVK

const ORDERS = Dict(parse(Int, a) => parse(Int, b)
for (a, b) in split.(split(get(ENV, "MORFE_BENCH_ORDERS", "2=>19,3=>9,4=>6,5=>5"), ","), "=>"))
const ROMS = parse.(Int, split(get(ENV, "MORFE_BENCH_ROMS", "2,3,4,5"), ","))
const LEVELS = 0:4
const PHASES = ("setup", "gradients", "stresses", "contraction", "assembly")

material = SVK.SVKMaterial(E = 160e3, ν = 0.22, ρ = 2.32e-3)
damping = SVK.RayleighDamping(α = 0.0, β = 0.0)
grid = generate_grid(Hexahedron, (60, 4, 2), Vec((0.0, 0.0, 0.0)), Vec((1000.0, 10.0, 24.0)))
addfacetset!(grid, "Dirichlet", x -> abs(x[1]) < 1e-8 || abs(x[1] - 1000) < 1e-8)
case = SVK.mechanical_model(grid; material, damping, dirichlet = "Dirichlet")
N = size(case.B[1], 1)

# ── M3: benchmark.jl's closure, stopping after phase LEVEL ───────────────────────────────
struct Closure{LEVEL, T, B, P, F}
    t::T
    ∇u::B
    P::P
    Fe::F
end
function Closure{LEVEL}(t) where {LEVEL}
    ∇u = similar(fem_qp_buffer(t), t.deg, fem_n_qp(t))
    P = Vector{typeof(fem_flux_zero(t))}(undef, fem_n_qp(t))
    Fe = zeros(ComplexF64, fem_ndofs_per_cell(t))
    return Closure{LEVEL, typeof(t), typeof(∇u), typeof(P), typeof(Fe)}(t, ∇u, P, Fe)
end

function (u::Closure{LEVEL})(res, args::Vararg{Any, N}) where {LEVEL, N}
    t = u.t
    n_qp = fem_n_qp(t)
    for element in fem_elements(t)
        fem_reinit!(element, t)
        if LEVEL >= 1
            for k in 1:N
                j = 0
                for i in 1:(k - 1)
                    args[i] === args[k] && (j = i; break)
                end
                j == 0 ? scatter_qp!(view(u.∇u, k, :), args[k], element, t) :
                copyto!(view(u.∇u, k, :), view(u.∇u, j, :))
            end
        end
        if LEVEL >= 2
            for q in 1:n_qp
                u.P[q] = flux_qp(ntuple(k -> u.∇u[k, q], Val(N)), element, q, t)
            end
        end
        if LEVEL >= 3
            fill!(u.Fe, 0)
            for q in 1:n_qp
                contract_qp!(u.Fe, u.P[q], element, q, fem_getdetJdV(element, q, t), t)
            end
        end
        LEVEL >= 4 && assemble_element!(res, u.Fe, element, t)
    end
    return res
end

# ── M4: MORFE's _replay_fem!, stopping after phase LEVEL ─────────────────────────────────
function single_pass!(result, W, cache, l, ::Val{LEVEL}) where {LEVEL}
    group = only(cache.fem_groups)   # the quadratic and cubic terms share one element loop
    cols = group.cols[l]
    (; terms, qp_buffer, Fe) = group
    external_arguments = cache.external_arguments
    entries = group.entries[l]
    driver = first(terms)
    flux0 = fem_flux_zero(driver)
    sink = flux0
    isempty(cols) && return sink
    n_qp = fem_n_qp(driver)
    for element in fem_elements(driver)
        fem_reinit!(element, driver)
        if LEVEL >= 1
            for (i, (order, col)) in enumerate(cols)
                scatter_qp!(@view(qp_buffer[i, 1:n_qp]), @view(W[:, order, col]), element, driver)
            end
        end
        LEVEL >= 3 && fill!(Fe, zero(eltype(Fe)))
        if LEVEL >= 2
            for q in 1:n_qp
                dΩ = fem_getdetJdV(element, q, driver)
                coefficient_values = @view qp_buffer[1:length(cols), q]
                flux = MORFE.accumulate_terms_qp!(flux0, Fe, coefficient_values,
                    external_arguments,
                    entries, terms, element, q, dΩ)
                LEVEL >= 3 ? contract_qp!(Fe, flux, element, q, dΩ, driver) : (sink += flux)
            end
        end
        LEVEL >= 4 && assemble_element!(result, Fe, element, driver)
    end
    return sink
end

partner(e, n) = Tuple(n >= 4 ? (e[2], e[1], e[4], e[3], e[5:end]...) : (e[2], e[1], e[3:end]...))
assembled(mset, n) = [l for l in eachindex(mset.exponents)
                      if sum(mset.exponents[l]) >= 2 &&
                         find_in_set(mset, partner(Tuple(mset.exponents[l]), n)) >= l]

terms = Tuple(case.term_factory(d) for d in (2, 3))
closure_model(level) = NthOrderModel(case.B,
    map(t -> MultilinearMap(Closure{level}(t); multiindex = t.multiindex,
            multiplicity_external = 0, fully_asymmetric = false), terms))
closures = Dict(level => closure_model(level) for level in LEVELS)
single = NthOrderModel(case.B, terms)

# The stresses of M4's level 2 are summed into a value kept here, so they cannot be dropped.
const SINK = Ref{Any}(nothing)

"""Time of one route at one level, summed over the monomials `ls`."""
function route_time(route, level, W, caches, ls, buf)
    t = 0.0
    sink = nothing
    for l in ls
        if route == "M3"
            t += @elapsed compute_multilinear_terms!(buf, closures[level], l, W, caches[level])
        else
            fill!(buf, 0)
            t += @elapsed(sink = single_pass!(buf, W.poly.coefficients, caches[:single], l,
                Val(level)))
        end
    end
    SINK[] = sink
    return t
end

caches_for(W) = merge(Dict{Any, Any}(level => build_multilinear_terms_cache(closures[level], W)
                                     for level in LEVELS),
    Dict{Any, Any}(:single => build_multilinear_terms_cache(single, W)))

# Compile every level of both routes on a small expansion first: the types do not depend on
# the multiindex set, so the timed sweeps run compiled code.
let mset = build_multiindex_set(3, 2), W = Parametrisation(
        DensePolynomial(randn(Xoshiro(0), ComplexF64, N, 2, length(mset)), mset), 0)
    caches = redirect_stderr(devnull) do
        caches_for(W)
    end
    buf = zeros(ComplexF64, N)
    for route in ("M3", "M4"), level in LEVELS
        route_time(route, level, W, caches, assembled(mset, 2), buf)
    end
end

out = joinpath(@__DIR__, "results", "data", "single_pass_phases.csv")
mkpath(dirname(out))
open(out, "w") do io
    println(io, "# $(now()), Julia $(VERSION), $(Sys.cpu_info()[1].model), ",
        "$(Threads.nthreads()) thread(s), medium beam mesh, $N dofs; time summed over the ",
        "assembled monomials; level k runs phases 0..k of ", join(PHASES, ", "))
    println(io, "n,order,monomials,route,level,time_s")
    for n in ROMS
        order = ORDERS[n]
        mset = build_multiindex_set(order, n)
        ls = assembled(mset, n)
        W = Parametrisation(DensePolynomial(zeros(ComplexF64, N, 2, length(mset)), mset), 0)
        randn!(Xoshiro(1), W.poly.coefficients)
        caches = redirect_stderr(devnull) do
            caches_for(W)
        end
        buf = zeros(ComplexF64, N)
        # the copied single pass, at full level, must reproduce MORFE's right-hand side
        reference = zeros(ComplexF64, N)
        compute_multilinear_terms!(reference, single, last(ls), W, caches[:single])
        fill!(buf, 0)
        single_pass!(buf, W.poly.coefficients, caches[:single], last(ls), Val(4))
        buf == reference || error("n = $n: the phase copy of the single pass deviates")
        GC.gc()
        times = Dict{Tuple{String, Int}, Float64}()
        for level in LEVELS, route in ("M3", "M4")   # the routes alternate at every level
            times[(route, level)] = route_time(route, level, W, caches, ls, buf)
            println(io, join((n, order, length(ls), route, level, times[(route, level)]), ","))
            flush(io)
        end
        for route in ("M3", "M4")
            phase = [times[(route, k)] - (k == 0 ? 0.0 : times[(route, k - 1)]) for k in LEVELS]
            @printf("n = %d, order %2d, %3d monomials, %s: total %8.2f s = %s\n", n, order,
                length(ls), route, times[(route, 4)],
                join((@sprintf("%s %.2f", p, t) for (p, t) in zip(PHASES, phase)), ", "))
        end
        W = caches = nothing
        GC.gc()
    end
end
println("wrote ", out)
