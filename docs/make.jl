using Documenter
using Graphics

makedocs(
    sitename = "Graphics",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [Graphics],
    pages = ["index.md", "reference.md"],
    clean = false,
    checkdocs = :exports
)

DocMeta.setdocmeta!(Graphics, :DocTestSetup, :(using Graphics); recursive=true)

deploydocs(
    repo = "github.com/JuliaGraphics/Graphics.jl.git"
)
