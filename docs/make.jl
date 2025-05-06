using Documenter
using DiscriminantVariety

makedocs(
    sitename = "DiscriminantVariety.jl",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    remotes = nothing
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#
deploydocs( 
    repo = "github.com/sumiya11/DiscriminantVariety.jl.git",
    push_preview = true
)

