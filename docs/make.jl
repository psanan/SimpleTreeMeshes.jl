push!(LOAD_PATH,"../src/")

using Documenter
using SimpleTreeMeshes


makedocs(;
         sitename="SimpleTreeMeshes.jl Documentation",
         modules  = [SimpleTreeMeshes],
         format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         authors = "Patrick Sanan",
         pages=[
                "Home" => "index.md"
               ],
)

deploydocs(;
         repo="github.com/psanan/SimpleTreeMeshes.jl.git",
         devbranch = "main",
         branch = "gh-pages",
         target = "build",
         push_preview = true,
)
