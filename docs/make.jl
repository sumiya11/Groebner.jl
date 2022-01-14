
push!(LOAD_PATH, "../src/")

using Documenter, FastGroebner

makedocs(
   sitename="FastGroebner",
   authors="sumiya11",
   modules=[FastGroebner],

   pages=[
      "Home" => "index.md",
   ]
)

deploydocs(
   branch = "gh-pages",
   repo = "github.com/sumiya11/FastGroebner.jl.git";
   push_preview = true
)
