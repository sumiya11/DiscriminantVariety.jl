using Documenter
using DiscriminantVariety

# Individual files
About = "Introduction" => "index.md"
Manual = "Manual" => [
    "Discriminant varieties" => "discriminant.md",
    "Computation over QQ" => "modular.md",
]
PublicAPI = "Public API" => "api.md"
Reference = "Reference" => [PublicAPI]

PAGES = [
    About,
    Manual,
    Reference
]

# Using MathJax3 for better math rendering
format = Documenter.HTML(
    mathengine = Documenter.MathJax3(),
    prettyurls = get(ENV, "CI", nothing) == "true"
)

makedocs(
    sitename = "DiscriminantVariety.jl",
    format = format,
    modules = [DiscriminantVariety],
    repo = "https://github.com/sumiya11/DiscriminantVariety.jl/blob/{commit}{path}#{line}",
    checkdocs = :exports,
    pages = PAGES
)

deploydocs(
    repo = "https://github.com/sumiya11/DiscriminantVariety.jl",
    push_preview = true
)
