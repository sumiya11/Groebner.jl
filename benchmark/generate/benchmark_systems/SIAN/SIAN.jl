import AbstractAlgebra

eqpath(path1, path2) = normpath(path1) == normpath(path2)

for (root, dir, files) in walkdir((@__DIR__))
    for file in files
        if eqpath((@__DIR__) * "/$file", (@__FILE__))
            continue
        end
        include((@__DIR__) * "/$file")
    end
end

function load_SIAN_all(; np=AbstractAlgebra, ground=np.QQ, internal_ordering=:degrevlex)
    systems = []
    for (root, dir, files) in walkdir((@__DIR__))
        for file in files
            if eqpath((@__DIR__) * "/$file", (@__FILE__))
                continue
            end
            name = Symbol(split(file, "/")[end][1:(end - 3)])
            sys = eval(:($name(np=($np), k=($ground), ordering=($(Meta.quot(ordering))))))
            push!(systems, (string(name), sys))
        end
    end
    systems
end
