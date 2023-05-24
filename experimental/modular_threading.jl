
using Base.Threads

function normalize_row(row)
    # @info "I am normalizing row" threadid()
    mult = inv(row[1])
    for i in eachindex(row)
        row[i] += row[i]
        row[i] *= mult
    end
    row
end

function uwu_spawn(rows)
    tasks = Vector{Task}(undef, 0)
    for row in rows
        task = @spawn normalize_row(row)
        push!(tasks, task)
    end
    for task in tasks
        wait(task)
    end
end

function uwu(rows)
    for row in rows
        normalize_row(row)
    end
end

rows = [randn(1000) for i in 1:20];
