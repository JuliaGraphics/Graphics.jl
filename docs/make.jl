using Documenter
using Graphics

makedocs(
    sitename = "Graphics",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [Graphics],
    pages = ["index.md", "reference.md"]
)

deploydocs(
    repo = "github.com/JuliaGraphics/Graphics.jl.git"
)
