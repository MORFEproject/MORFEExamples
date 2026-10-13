# Draws results/data/batching_benchmark.csv, written by benchmark.jl, against the number of
# solved monomials, with dashed lines at the degree boundaries:
#   1. cumulative right-hand-side time and allocated memory, closure against batched;
#   2. cumulative mesh traversals.
#   julia plot_benchmark.jl
# Writes results/figures/batching_benchmark.png and, for the MORFE.jl paper's pgfplots
# figure, results/data/batching_<mesh>.dat (cumulative series, one per mesh) and
# results/data/batching_degrees.tex (the degree boundaries as TeX macros).
using Pkg: Pkg
Pkg.activate(dirname(@__DIR__))
using CairoMakie, Printf

data = joinpath(@__DIR__, "results", "data")
lines_ = readlines(joinpath(data, "batching_benchmark.csv"))
comment = filter(startswith("#"), lines_)
body = filter(!startswith("#"), lines_)
header = Symbol.(split(body[1], ','))
rows = [NamedTuple{Tuple(header)}(Tuple(split(l, ','))) for l in body[2:end]]

meshes = unique(r.mesh for r in rows)
series(mesh, route) = let rs = sort(filter(r -> r.mesh == mesh && r.route == route, rows);
        by = r -> parse(Int, r.k))
    (k = parse.(Int, getfield.(rs, :k)), degree = parse.(Int, getfield.(rs, :degree)),
        passes = cumsum(parse.(Int, getfield.(rs, :passes))),
        time = cumsum(parse.(Float64, getfield.(rs, :rhs_time_s))),
        mb = cumsum(parse.(Int, getfield.(rs, :rhs_alloc_bytes))) ./ 2^20,
        dofs = parse(Int, rs[1].n_dofs))
end
S = Dict((m, r) => series(m, r) for m in meshes, r in ("closure", "batched"))

# Degree boundaries, from the solve order (the same for every mesh).
ref = S[(meshes[1], "closure")]
cuts = [ref.k[i] + 0.5 for i in 1:(length(ref.k) - 1) if ref.degree[i + 1] != ref.degree[i]]
spans = [(d, (minimum(ref.k[ref.degree .== d]) + maximum(ref.k[ref.degree .== d])) / 2)
         for d in unique(ref.degree)]

# ── quick-look figure ────────────────────────────────────────────────────────────────
shade(base, i) = RGBAf(Makie.to_color(base).r, Makie.to_color(base).g,
    Makie.to_color(base).b, 0.35 + 0.65 * i / length(meshes))
fig = Figure(size = (1200, 380))
labels = ("right-hand-side time [s]", "allocated memory [MiB]", "mesh traversals")
# round 1-2-5 ticks on the log axes rather than Makie's 10^x.5 labels
ticks125(lo, hi) = filter(v -> lo <= v <= hi,
    [m * 10.0^e for e in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in (1, 2, 5)])
logticks(vals) = (t = ticks125(minimum(vals), maximum(vals)); (t, [@sprintf("%g", v) for v in t]))
all_time = [v for s in values(S) for v in s.time]
all_passes = [v for s in values(S) for v in s.passes]
yt = (logticks(all_time), Makie.automatic, logticks(all_passes))
axes = [Axis(fig[1, j]; xlabel = "solved monomials", ylabel = labels[j],
            yscale = j == 2 ? identity : log10, yticks = yt[j], xticks = 0:5:ref.k[end])
        for j in 1:3]
for ax in axes
    vlines!(ax, cuts; color = (:gray, 0.6), linestyle = :dash, linewidth = 1)
end
ytop = 2 * maximum(S[(m, "closure")].time[end] for m in meshes)   # room for the labels
for (d, mid) in spans
    text!(axes[1], mid, ytop; text = "degree $d", align = (:center, :top), fontsize = 10,
        color = :gray, rotation = 0)
end
for (i, m) in enumerate(meshes), (route, base) in (("closure", :darkorange), ("batched", :royalblue))
    s = S[(m, route)]
    label = "$route, $m ($(s.dofs) dofs)"
    lines!(axes[1], s.k, s.time; color = shade(base, i), label)
    lines!(axes[2], s.k, s.mb; color = shade(base, i))
    i == 1 && lines!(axes[3], s.k, s.passes; color = base, label = route)
end
ylims!(axes[1], nothing, 1.2 * ytop)
Legend(fig[1, 4], axes[1]; framevisible = false, labelsize = 11)
axislegend(axes[3]; position = :lt, framevisible = false)
isempty(comment) || Label(fig[0, 1:4], first(comment); fontsize = 10, color = :gray)
mkpath(joinpath(@__DIR__, "results", "figures"))
save(joinpath(@__DIR__, "results", "figures", "batching_benchmark.png"), fig; px_per_unit = 2)

# ── pgfplots exports ─────────────────────────────────────────────────────────────────
for m in meshes
    c, b = S[(m, "closure")], S[(m, "batched")]
    open(joinpath(data, "batching_$(m).dat"), "w") do io
        println(io, "# cumulative over the solved monomials; ", m, ", ", c.dofs, " dofs")
        println(io, "k degree passes_closure passes_batched time_closure time_batched ",
            "mib_closure mib_batched")
        for i in eachindex(c.k)
            @printf(io, "%d %d %d %d %.6g %.6g %.6g %.6g\n", c.k[i], c.degree[i],
                c.passes[i], b.passes[i], c.time[i], b.time[i], c.mb[i], b.mb[i])
        end
    end
end
open(joinpath(data, "batching_degrees.tex"), "w") do io
    # Literal TikZ commands rather than loops: pgfplots defers what is drawn inside an axis,
    # and loop variables would be gone by then.
    println(io, "% Degree boundaries of the solve order, written by plot_benchmark.jl.")
    println(io, "\\def\\BatchingDegreeCuts{", join(("\\draw[degreecut] " *
                "({axis cs:$x,1} |- {rel axis cs:0,0}) -- ({axis cs:$x,1} |- {rel axis cs:0,1});"
                for x in cuts)), "}")
    println(io, "\\def\\BatchingDegreeLabels{", join(("\\node[degreelabel] at " *
                "({axis cs:$mid,1} |- {rel axis cs:0,1}) {$d};" for (d, mid) in spans)), "}")
    println(io, "\\def\\BatchingLastMonomial{", ref.k[end], "}")
end
println("wrote batching_benchmark.png, ", join(("batching_$(m).dat" for m in meshes), ", "),
    " and batching_degrees.tex")
