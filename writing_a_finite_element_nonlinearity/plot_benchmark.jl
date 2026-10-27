# Draws results/data/single_pass_benchmark.csv, written by benchmark.jl, against the number of
# monomials of the expansion, in the solver's order. The solver assembles a right-hand side
# for one monomial of each complex-conjugate pair and fills the other by conjugation, and
# the linear monomials need none, so the cumulative curves stay flat across those.
#   1. speed-up of the single-pass route, closure time / single-pass time;
#   2. the bytes the closure allocates during assembly; the single pass allocates none,
#      and its total, if any, is in the axis label;
#   3. the expansion order reached, the same for every mesh: order d at the last monomial of
#      degree d, from order 1, joined by straight lines in log-log axes;
#   4. the factorisations evaluated, the same for every mesh (one closure call each).
# Colour: reduced dimension n; solid lines for the medium mesh, and dashed and dotted for any
# other mesh in the CSV. Markers end each degree and cycle by order (circle, star, triangle, x,
# square) for easy counting; in panels 1, 2 and 4, a faint grey line joins the markers of
# the same order across n. The panels sit in a 2 x 2 grid. The monomial axis is
# logarithmic, labelled at the decades; so are the allocation, order and factorisation axes.
#   julia plot_benchmark.jl
# Writes results/figures/single_pass_benchmark.png and, for the MORFE.jl paper's pgfplots
# figure (names as in its figures/data/), results/data/single_pass_assembly_<mesh>_n<n>.dat,
# one cumulative series per mesh and n, results/data/single_pass_assembly_orders_n<n>.dat,
# the degree boundaries per n, and results/data/single_pass_assembly_<mesh>_isoorders.dat,
# the same-order markers across n.
using Pkg: Pkg
Pkg.activate(dirname(@__DIR__))
using CairoMakie, Printf

data = joinpath(@__DIR__, "results", "data")
lines_ = readlines(joinpath(data, "single_pass_benchmark.csv"))
comment = filter(startswith("#"), lines_)
body = filter(!startswith("#"), lines_)
header = Symbol.(split(body[1], ','))
rows = [NamedTuple{Tuple(header)}(Tuple(split(l, ','))) for l in body[2:end]]

meshes = unique(r.mesh for r in rows)
roms = sort(unique(parse(Int, r.n) for r in rows))
num(rs, f, T = Float64) = parse.(T, getfield.(rs, f))

# Position of the last monomial of each degree 1:order in the graded set of n variables.
degree_ends(n, order) = cumsum(binomial(n + d - 1, d) for d in 1:order)

function series(mesh, n)
    pick(route) = sort(filter(r -> r.mesh == mesh && parse(Int, r.n) == n && r.route == route, rows);
        by = r -> parse(Int, r.k))
    c, b = pick("closure"), pick("single_pass")
    order = parse(Int, c[1].order)
    ends = degree_ends(n, order)
    idx = num(c, :monomial_idx, Int)
    @assert idx == num(b, :monomial_idx, Int)
    # Per-monomial cost over the whole set, zero where nothing is assembled.
    spread(rs, f) = (v = zeros(ends[end]); v[idx] = num(rs, f); cumsum(v))
    Tc, Tb = spread(c, :rhs_time_s), spread(b, :rhs_time_s)
    Bc, Bb = spread(c, :rhs_alloc_bytes), spread(b, :rhs_alloc_bytes)
    F = spread(c, :passes)   # the closure evaluates each factorisation once
    ls = first(idx):ends[end]
    return (k = collect(ls), degree = [findfirst(>=(l), ends) for l in ls],
        last = [l in ends for l in ls], speedup = Tc[ls] ./ Tb[ls],
        allocated = Bc[ls], allocated_single = Bb[end], factorisations = F[ls],
        ends, order,
        dofs = parse(Int, c[1].n_dofs))
end
S = Dict((m, n) => series(m, n) for m in meshes for n in roms)

# The order reached: order d at the last monomial of degree d, from order 1 at the last
# linear monomial. Joining them gives a constant slope within each order in log-log axes.
boundaries(s) = (k = s.ends, order = collect(1:s.order))

# Markers cycle by order: circle, star, triangle, x, square, then again from order 6.
markers = (:circle, :star5, :utriangle, :xcross, :rect)
sizes = (8, 13, 10, 12, 7)   # balanced by eye: star and cross run small, the square large
marker(d) = markers[mod(d - 1, 5) + 1]
markersize(d) = sizes[mod(d - 1, 5) + 1]

# Same-order markers across n, per mesh: for every degree, the points of the curves that
# reach it, by increasing n.
function isoorders(m)
    top = maximum(S[(m, n)].order for n in roms)
    return [(d, [(S[(m, n)].ends[d], i) for n in roms if 2 <= d <= S[(m, n)].order
                 for i in findall(==(S[(m, n)].ends[d]), S[(m, n)].k)]) for d in 2:top]
end

# ── quick-look figure ────────────────────────────────────────────────────────────────
colours = Dict(2 => "#1f77c8", 3 => "#14a56e", 4 => "#f58214", 5 => "#b4325a")
styles = Dict("medium" => :solid, "large" => :dash, "extra-large" => :dot)
fig = Figure(size = (1150, 820))
single_bytes = maximum(s.allocated_single for s in values(S))
labels = ("speed-up  (closure / single-pass time)",
    "bytes allocated by the closure  (single pass: $(round(Int, single_bytes)))",
    "order", "factorisations")
fields = (:speedup, :allocated)
top = maximum(S[(meshes[1], n)].order for n in roms)
# Logarithmic monomial axis: labels at the decades, tick marks at every 1-9 multiple.
logx = (xscale = log10, limits = ((1, 300), nothing),
    xticks = ([1, 10, 100], ["1", "10", "100"]), xminorticksvisible = true,
    xminorticks = [2:9; 20:10:90; 200:100:300], xgridvisible = false)
axes = [Axis(fig[1, 1]; ylabel = labels[1], ygridvisible = false, logx...),
    Axis(fig[1, 2]; ylabel = labels[2], ygridvisible = false,
        yscale = log10,
        yticks = ([1e3, 1e4, 1e5, 1e6, 1e7], ["10³", "10⁴", "10⁵", "10⁶", "10⁷"]),
        yminorticksvisible = true,
        yminorticks = [k * 10.0^e for e in 3:6 for k in 2:9],
        logx...),
    Axis(fig[2, 1]; xlabel = "monomials", ylabel = labels[3], yscale = log10,
        yticks = [1, 2, 3, 5, 10, 19], yminorticks = 1:top, yminorgridvisible = true,
        yminorgridcolor = (:black, 0.12), logx...),
    Axis(fig[2, 2]; xlabel = "monomials", ylabel = labels[4], yscale = log10,
        yticks = ([1, 10, 100, 1000, 10000], ["1", "10", "10²", "10³", "10⁴"]),
        logx...)]   # horizontal gridlines at each decade
for m in meshes, (d, pts) in isoorders(m), (j, f) in ((1, :speedup), (2, :allocated), (4, :factorisations))
    m == meshes[1] || j != 4 || continue   # factorisations do not depend on the mesh
    lines!(axes[j], first.(pts), [getfield(S[(m, n)], f)[i]
                                  for (n, (_, i)) in zip(filter(n -> 2 <= d <= S[(m, n)].order, roms), pts)];
        color = (:gray, 0.2), linewidth = 1)
end
for n in roms, m in meshes, j in 1:2
    s = S[(m, n)]
    lines!(axes[j], s.k, getfield(s, fields[j]); color = colours[n], linestyle = styles[m],
        linewidth = 2.5)
    for i in findall(s.last)
        scatter!(axes[j], [s.k[i]], [getfield(s, fields[j])[i]]; color = colours[n],
            marker = marker(s.degree[i]), markersize = markersize(s.degree[i]))
    end
end
for n in roms
    s = S[(meshes[1], n)]
    lines!(axes[4], s.k, s.factorisations; color = colours[n], linewidth = 2.5)
    for i in findall(s.last)
        scatter!(axes[4], [s.k[i]], [s.factorisations[i]]; color = colours[n],
            marker = marker(s.degree[i]), markersize = markersize(s.degree[i]))
    end
end
for n in roms
    b = boundaries(S[(meshes[1], n)])
    lines!(axes[3], b.k, b.order; color = colours[n], linewidth = 2.5)
    for i in 1:length(b.k)
        scatter!(axes[3], [b.k[i]], [b.order[i]]; color = colours[n], marker = marker(b.order[i]),
            markersize = markersize(b.order[i]))
    end
end
meshlegend = length(meshes) > 1   # line styles only matter when there are several meshes
Legend(fig[1:2, 3],
    [[LineElement(color = colours[n]) for n in roms];
        meshlegend ? [LineElement(color = :gray, linestyle = styles[m]) for m in meshes] : []],
    [["n = $n, order $(S[(meshes[1], n)].order)" for n in roms];
        meshlegend ? ["$m ($(S[(m, roms[1])].dofs) dofs)" for m in meshes] : String[]];
    framevisible = false, labelsize = 11)
rowgap!(fig.layout, 1, 8)   # the top row shares the monomial axis of the bottom row
isempty(comment) || Label(fig[0, 1:3], first(comment); fontsize = 10, color = :gray)
mkpath(joinpath(@__DIR__, "results", "figures"))
save(joinpath(@__DIR__, "results", "figures", "single_pass_benchmark.png"), fig; px_per_unit = 2)

# ── pgfplots exports ─────────────────────────────────────────────────────────────────
for ((m, n), s) in S
    open(joinpath(data, "single_pass_assembly_$(m)_n$(n).dat"), "w") do io
        println(io, "# cumulative over the monomials of the expansion, in solve order; ",
            "allocated: bytes the closure allocates; ", m,
            ", ", s.dofs, " dofs, n = ", n, ", order ", s.order, "; last = 1 ends a degree")
        println(io, "k degree last speedup allocated factorisations")
        for i in eachindex(s.k)
            @printf(io, "%d %d %d %.6g %.6g %d\n", s.k[i], s.degree[i], s.last[i],
                s.speedup[i], s.allocated[i], s.factorisations[i])
        end
    end
end
for n in roms
    b = boundaries(S[(meshes[1], n)])
    open(joinpath(data, "single_pass_assembly_orders_n$(n).dat"), "w") do io
        println(io, "# degree boundaries, n = ", n, ": order d at the last monomial of ",
            "degree d")
        println(io, "k order")
        foreach(i -> println(io, b.k[i], " ", b.order[i]), eachindex(b.k))
    end
end
for m in meshes
    open(joinpath(data, "single_pass_assembly_$(m)_isoorders.dat"), "w") do io
        println(io, "# same-order markers across n, one block per degree (blank line between); ", m)
        println(io, "degree n k speedup allocated factorisations")
        for (d, pts) in isoorders(m)
            for (n, (k, i)) in zip(filter(n -> 2 <= d <= S[(m, n)].order, roms), pts)
                s = S[(m, n)]
                @printf(io, "%d %d %d %.6g %.6g %d\n", d, n, k, s.speedup[i], s.allocated[i],
                    s.factorisations[i])
            end
            println(io)
        end
    end
end
println("wrote single_pass_benchmark.png, ", length(S), " single_pass_assembly_<mesh>_n<n>.dat, ",
    length(roms), " single_pass_assembly_orders_n<n>.dat and ", length(meshes),
    " single_pass_assembly_<mesh>_isoorders.dat files")
