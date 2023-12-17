include("stuff2.lj")

struct Tre
    x::Union{Tre, Nothing}
end

y = Tre(1)
