using Groebner
using Documenter

makedocs(
    modules = [Groebner],
    sitename = "Groebner.jl",
    doctest = true,
    linkcheck = true,
    checkdocs = :exports,
    warnonly=true,
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Interface" => "interface.md",
    ],
)


deploydocs(
    repo   = "github.com/sumiya11/Groebner.jl.git"
)
