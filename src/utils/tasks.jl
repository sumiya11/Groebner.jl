# This file is a part of Groebner.jl. License is GNU GPL v2.

function split_round_robin(data, tasks::Int)
    chunk_size = max(1, div(length(data), tasks, RoundUp))
    data_chunks = [
        [i + tasks * (j - 1) for j in 1:chunk_size if i + tasks * (j - 1) <= length(data)] for
        i in 1:tasks
    ]
    data_chunks
end

@assert split_round_robin(1:10, 3) == [[1, 4, 7, 10], [2, 5, 8], [3, 6, 9]]
