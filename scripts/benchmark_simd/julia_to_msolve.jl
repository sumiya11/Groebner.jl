using Groebner, Nemo

function to_msolve(sys)
  s = ""
  s = s*join(map(string, gens(parent(sys[1]))), ",")*"\n"
  s = s*string(characteristic(parent(sys[1])))*"\n"
  s = s*join(map(string, sys), ",\n")
  s
end

systems = [
  ("cyclic7_zp.txt", Groebner.Examples.cyclicn(7,k=GF(2^30+3))),
  ("cyclic7_qq.txt", Groebner.Examples.cyclicn(7))
]

for (name, sys) in systems
  sys = to_msolve(sys)
  io = open(joinpath(@__DIR__, name), "w")
  println(io, sys)
  close(io)
end
