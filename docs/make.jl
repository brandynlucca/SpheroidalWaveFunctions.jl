using Documenter
using SpheroidalWaves

makedocs(
    sitename = "SpheroidalWaves.jl",
    modules = [SpheroidalWaves],
    authors = "SpheroidalWaves contributors",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Math and Usage" => "math-and-usage.md",
        "Backend Overrides" => "backend-overrides.md",
    ],
)

deploydocs(
    repo = "github.com/brandynlucca/SpheroidalWaves.jl.git",
    devbranch = "main",
)

