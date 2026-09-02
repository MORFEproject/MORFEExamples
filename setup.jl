using Pkg

Pkg.activate(@__DIR__)
Pkg.add("MORFE")
Pkg.add(url="https://github.com/MORFEproject/MORFEFerrite.jl")
Pkg.add("IJulia")

Pkg.instantiate()

using IJulia
installkernel(
    "MORFEExamples",
    env=Dict("JULIA_PROJECT" => @__DIR__)
)