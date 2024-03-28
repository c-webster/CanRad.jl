push!(LOAD_PATH,"../src/")
using CanRad
using Documenter
makedocs(
         sitename = "CanRad.jl",
         modules  = [CanRad],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/c-webster/CanRad.jl",
)