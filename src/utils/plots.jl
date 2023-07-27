# Partially adapted from and entirely inspired by UnicodePlots.jl 

const UnicodeType = UInt32

mutable struct Canvas
    grid::Matrix{UnicodeType}
    width::Int
    height::Int
    x_origin::Int
    y_origin::Int
    x_scale::Any
    y_scale::Any

    function Canvas(xlims, ylims; width=10, height=10)
        xmin, xmax = xlims
        ymin, ymax = ylims
        grid = fill(Char(0x0020), width, height)
        x_scale = width / (xmax - xmin + 1)
        y_scale = height / (ymax - ymin + 1)
        x_origin = 1
        y_origin = 1
        new(grid, width, height, x_origin, y_origin, x_scale, y_scale)
    end
end

function point!(c::Canvas, x::Integer, y::Integer)
    x_pixel, y_pixel = to_pixel(c, x, y)
    pixel!(c, x_pixel, y_pixel)
end

function to_pixel(c::Canvas, x, y)
    x_pixel = max(1, round(Int, x * c.x_scale - c.x_scale))
    y_pixel = max(1, round(Int, y * c.y_scale - c.x_scale))
    x_pixel, y_pixel
end

function pixel_to_char_point(c::Canvas, x_pixel, y_pixel)
    x_char = x_pixel
    y_char = y_pixel
    x_char, y_char
end

function pixel!(c::Canvas, x_pixel::Integer, y_pixel::Integer)
    char_x, char_y = pixel_to_char_point(c, x_pixel, y_pixel)
    offset_x, offset_y = floor(Int, c.x_scale), floor(Int, c.y_scale)
    for i in char_x:(char_x + offset_x)
        for j in char_y:(char_y + offset_y)
            val = c.grid[i, j]
            c.grid[i, j] = UInt('.')
        end
    end
    true
end

Base.show(io::IO, c::Canvas) = _show(io, c)

function _show(io::IO, c::Canvas)
    for i in 1:(c.width)
        for j in 1:(c.height)
            print(io, Char(c.grid[i, j]))
        end
        println(io)
    end
    nothing
end