
push!(LOAD_PATH, "../src/")

using Documenter, Groebner

makedocs(
   sitename="Groebner.jl",
   authors="sumiya11",
   modules=[Groebner],

   pages=[
      "Home" => "index.md",
      "Tutorials" => Any[
            "tutorials/system_solving.md",
       ],
       "Notes on theory" => "theory.md"
   ]
)

deploydocs(
   repo = "github.com/sumiya11/Groebner.jl.git",
   branch = "main",
   devbranch = "main",
   push_preview = true
)
