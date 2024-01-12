import AbstractAlgebra

for (root, dir, files) in walkdir((@__DIR__))
    for file in files
        if (@__DIR__) * "/$file" == (@__FILE__)
            continue
        end
        include((@__DIR__) * "/$file")
    end
end

function load_SIAN_all(; np=AbstractAlgebra, ground=np.QQ, ordering=:degrevlex)
    systems = []
    for (root, dir, files) in walkdir((@__DIR__))
        for file in files
            if (@__DIR__) * "/$file" == (@__FILE__)
                continue
            end
            name = Symbol(split(file, "/")[end][1:(end - 3)])
            sys = eval(:($name(np=$np, k=$ground, ordering=$(Meta.quot(ordering)))))
            push!(systems, (string(name), sys))
        end
    end
    systems
end
