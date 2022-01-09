
push!(LOAD_PATH, "../src/")

using Documenter, GroebnerBases

makedocs(
   sitename="GroebnerBases",
   authors="sumiya11",
   modules=[GroebnerBases],

   pages=[
      "Home" => "index.md",
   ]
)

deploydocs(
   repo = "github.com/sumiya11/GroebnerBases.git";
   push_preview = true,
   devbranch = "main"
)
