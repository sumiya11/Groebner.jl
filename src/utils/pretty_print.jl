# This file is a part of Groebner.jl. License is GNU GPL v2.

# Adapted from https://code.tecosaur.net/tec/About.jl.
# The license is MIT. Copyright (c) 2024 TEC.
function columnlist(
    io::IO,
    entries::Vector{<:AbstractString};
    maxcols::Int=8,
    maxwidth::Int=last(displaysize(io)),
    prefix::AbstractString="• ",
    spacing::Int=2
)
    isempty(entries) && return
    thecolumns = Vector{eltype(entries)}[]
    thecolwidths = Int[]
    for ncols in 2:maxcols
        columns = Vector{eltype(entries)}[]
        for col in Iterators.partition(entries, div(length(entries), ncols, RoundUp))
            push!(columns, collect(col))
        end
        widths = map(col -> map(textwidth, col), columns)
        colwidths = map(maximum, widths)
        layoutwidth = sum(colwidths) + ncols * textwidth(prefix) + (ncols - 1) * spacing
        if layoutwidth > maxwidth
            columns = Vector{eltype(entries)}[]
            for col in
                Iterators.partition(entries, div(length(entries), ncols - 1, RoundUp))
                push!(columns, collect(col))
            end
            widths = map(col -> map(textwidth, col), columns)
            colwidths = map(maximum, widths)
            thecolumns, thecolwidths = columns, colwidths
            break
        else
            thecolumns, thecolwidths = columns, colwidths
        end
    end
    for rnum in 1:length(first(thecolumns))
        for cnum in 1:length(thecolumns)
            rnum > length(thecolumns[cnum]) && continue
            cnum > 1 && print(io, ' '^spacing)
            print(io, prefix, rpad(thecolumns[cnum][rnum], thecolwidths[cnum]))
        end
        println(io)
    end
end
