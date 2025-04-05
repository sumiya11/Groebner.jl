# Adapted from JUMP. The license is MPL version 2.0.

using Groebner
using Documenter

# Setup

const _IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"

latex_platform = _IS_GITHUB_ACTIONS ? "docker" : "native"

# Documentation structure

PAGES = [
    "Home"      => "index.md", 
    "Examples"  => "examples.md", 
    "Interface" => "interface.md"
]

# Build the HTML documentation

# Needed to make Documenter think that there is a PDF in the right place when
# link checking. Inn production we replace this by running the LaTeX build.
write(joinpath(@__DIR__, "src", "Groebner.jl.pdf"), "")

@time Documenter.makedocs(
    modules   = [Groebner],
    sitename  = "Groebner.jl",
    doctest   = true,
    linkcheck = true,
    checkdocs = :exports,
    warnonly  = true,
    pages     = PAGES
)

# Build the PDF documentation

@time Documenter.makedocs(
    sitename = "Groebner.jl",
    authors  = "The Groebner.jl developers and contributors",
    format   = Documenter.LaTeX(; platform=latex_platform),
    build    = "latex_build",
    pages    = PAGES,
    debug    = true
)
cp(
    joinpath(@__DIR__, "latex_build", "Groebner.jl.pdf"),
    joinpath(@__DIR__, "build", "Groebner.jl.pdf");
    force=true
)

# Deploy build/

deploydocs(repo="github.com/sumiya11/Groebner.jl.git")
