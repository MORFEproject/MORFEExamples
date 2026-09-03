using Pkg

Pkg.activate(@__DIR__)
Pkg.add(url="https://github.com/MORFEproject/MORFEFerrite.jl")
Pkg.add("MORFE")
Pkg.add("IJulia")
Pkg.add("DifferentialEquations")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("Dates")
Pkg.add("Symbolics")
Pkg.add("StaticArrays")

Pkg.instantiate()

using IJulia
installkernel(
    "MORFEExamples",
    env=Dict("JULIA_PROJECT" => @__DIR__)
)