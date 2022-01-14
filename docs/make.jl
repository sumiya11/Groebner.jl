
push!(LOAD_PATH, "../src/")

using Documenter, FastGroebner

makedocs(
   sitename="FastGroebner.jl",
   authors="sumiya11",
   modules=[FastGroebner],

   pages=[
      "Home" => "index.md",
   ]
)

deploydocs(
   repo = "github.com/sumiya11/FastGroebner.jl.git",
   branch = "main",
   devbranch = "main",
   push_preview = true
)
