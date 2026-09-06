using Pkg

Pkg.activate(@__DIR__)
# Both packages come from their development branches rather than the registry.
# MORFEFerrite is not registered at all, and the registered MORFE v0.1.0 predates
# the post-processing API these examples use (observable_polynomial,
# cycle_amplitude). Switch MORFE back to Pkg.add("MORFE") once a release carries it.
Pkg.add(url="https://github.com/MORFEproject/MORFEFerrite.jl")
Pkg.add(url="https://github.com/MORFEproject/MORFE.jl")
Pkg.add("IJulia")
Pkg.add("DifferentialEquations")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("CairoMakie")
Pkg.add("Dates")
Pkg.add("Symbolics")
Pkg.add("StaticArrays")

Pkg.instantiate()

using IJulia
installkernel(
    "MORFEExamples",
    env=Dict("JULIA_PROJECT" => @__DIR__)
)