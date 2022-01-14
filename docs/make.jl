
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
   repo = "github.com/sumiya11/FastGroebner.git";
   push_preview = true,
   devbranch = "main"
)
