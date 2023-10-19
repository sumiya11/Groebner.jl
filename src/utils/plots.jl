# Automatic plotting of the sparsity pattern of F4 matrices

# Adapted from and entirely inspired by UnicodePlots.jl.
# The licence is MIT.
# https://github.com/JuliaPlots/UnicodePlots.jl/blob/master/LICENSE.md

const UnicodeType = UInt32

const BLANK_BRAILLE = 0x2800
const FULL_BRAILLE = 0x28ff

# braille dots composing ⣿
const BRAILLE_SIGNS = UnicodeType.([
    '⠁' '⠈'
    '⠂' '⠐'
    '⠄' '⠠'
    '⡀' '⢀'
])

#! format: off
const BORDER_SOLID = (
    tl = '┌',
    tr = '┐',
    bl = '└',
    br = '┘',
    t  = '─',
    l  = '│',
    b  = '─',
    r  = '│',
    c  = '+',
)
#! format: on

const ASPECT_RATIO = Ref(4 / 3)

# default display size for the default Canvas 
# (which has aspect ratio = 2) ==> (40, 15)
const DEFAULT_HEIGHT = Ref(15)
const DEFAULT_WIDTH = Ref(round(Int, DEFAULT_HEIGHT[] * 2ASPECT_RATIO[]))

y_pixel_per_char() = 4
x_pixel_per_char() = 2

out_stream_size(out_stream::Nothing) = displaysize()
out_stream_size(out_stream::IO) = displaysize(out_stream)

mutable struct Canvas
    grid::Matrix{UnicodeType}
    nrows::Int
    ncols::Int
    pixel_width::Int
    pixel_height::Int
    width::Float64
    height::Float64
    origin_x::Int
    origin_y::Int

    function Canvas(
        nrows::Integer,
        ncols::Integer;
        out_stream=stdout,
        height::Union{Integer, Nothing}=nothing,
        width::Union{Integer, Nothing}=nothing,
        max_height::Union{Integer, Nothing}=nothing,
        max_width::Union{Integer, Nothing}=nothing,
        fix_ar::Bool=false
    )
        char_height, char_width = get_canvas_dimensions(
            nrows,
            ncols,
            height,
            width,
            max_height,
            max_width,
            out_stream,
            fix_ar
        )
        char_height = max(char_height, 2)
        char_width = max(char_width, 5)
        pixel_height = char_height * y_pixel_per_char()
        pixel_width = char_width * x_pixel_per_char()
        grid = fill(BLANK_BRAILLE, char_height, char_width)
        origin_y = 0.0
        origin_x = 0.0
        height = 1.0 + nrows
        width = 1.0 + ncols
        new(
            grid,
            nrows,
            ncols,
            pixel_width,
            pixel_height,
            width,
            height,
            origin_x,
            origin_y
        )
    end
end

function get_canvas_dimensions(
    nrows::Integer,
    ncols::Integer,
    height::Union{Integer, Nothing},
    width::Union{Integer, Nothing},
    max_height::Union{Integer, Nothing},
    max_width::Union{Integer, Nothing},
    out_stream::IO,
    fix_ar::Bool
)
    canv_height = nrows / y_pixel_per_char()
    canv_width  = ncols / x_pixel_per_char()
    canv_ar     = canv_width / canv_height

    # min_canv_height := minimal number of y canvas characters
    # (holding y_pixel_per_char pixels) to represent the input data
    min_canv_height = ceil(Int, canv_height)
    min_canv_width  = ceil(Int, canv_width)

    height_diff = 0
    width_diff = 0

    term_height, term_width = out_stream_size(out_stream)
    max_height = max_height !== nothing ? max_height : term_height - height_diff
    max_width = max_width !== nothing ? max_width : term_width - width_diff

    (nrows == 0 && ncols == 0) && return 0, 0, max_width, max_height

    if width === nothing && height === nothing
        if min_canv_height > min_canv_width
            # long matrix (according to pixel density)
            width  = min(min_canv_height * canv_ar, max_width)
            height = min(width / canv_ar, max_height)
            width  = min(height * canv_ar, max_width)
        else
            # wide matrix
            height = min(min_canv_width / canv_ar, max_height)
            width  = min(height * canv_ar, max_width)
            height = min(width / canv_ar, max_height)
        end
    end

    if width ≡ nothing && height > 0
        width = min(height * canv_ar, max_width)
    elseif height ≡ nothing && width > 0
        height = min(width / canv_ar, max_height)
    end

    height = round(Int, height / (fix_ar ? ASPECT_RATIO[] : 1))  # optional terminal aspect ratio (4:3) correction
    width  = round(Int, width)

    # the canvas will target a (height, width) grid to represent the input data
    height, width, max_height, max_width
end

function point!(c::Canvas, x::Integer, y::Integer)
    x_pixel, y_pixel = point_to_pixel(c, x, y)
    pixel!(c, x_pixel, y_pixel)
end

scale_y_to_pixel(c::Canvas, y::Number) = (1 - (y - c.origin_y) / c.height) * c.pixel_height
scale_x_to_pixel(c::Canvas, x::Number) = (x - c.origin_x) / c.width * c.pixel_width

valid_y_pixel(c::Canvas, pixel_y::Integer) = 0 <= pixel_y <= c.pixel_height
valid_x_pixel(c::Canvas, pixel_x::Integer) = 0 <= pixel_x <= c.pixel_width

function point_to_pixel(c::Canvas, x::Integer, y::Integer)
    # NOTE the reverse of coordinates!
    x, y = y, x
    y = c.nrows - y + 1
    x_pixel = floor(Int, scale_x_to_pixel(c, x))
    y_pixel = floor(Int, scale_y_to_pixel(c, y))
    x_pixel, y_pixel
end

function pixel_to_char_point_off(c::Canvas, pixel_x::Integer, pixel_y::Integer)
    pixel_x >= c.pixel_width && (pixel_x -= 1)
    pixel_y >= c.pixel_height && (pixel_y -= 1)
    qx, rx = divrem(pixel_x, x_pixel_per_char())
    qy, ry = divrem(pixel_y, y_pixel_per_char())
    (Int(qx) + 1, Int(qy) + 1, Int(rx) + 1, Int(ry) + 1)
end

function pixel!(c::Canvas, pixel_x::Integer, pixel_y::Integer)
    valid_x_pixel(c, pixel_x) || return c
    valid_y_pixel(c, pixel_y) || return c
    char_x, char_y, char_x_off, char_y_off = pixel_to_char_point_off(c, pixel_x, pixel_y)
    if checkbounds(Bool, c.grid, char_y, char_x)
        if BLANK_BRAILLE ≤ (val = c.grid[char_y, char_x]) ≤ FULL_BRAILLE
            c.grid[char_y, char_x] = val | BRAILLE_SIGNS[char_y_off, char_x_off]
        end
    end
    c
end

Base.show(io::IO, c::Canvas) = _show(io, c)

function printrow(io::IO, row::AbstractVector)
    for j in 1:size(row, 1)
        print(io, Char(row[j]))
    end
end

function _show(io::IO, c::Canvas)
    for i in 1:size(c.grid, 1)
        printrow(io, c.grid[i, :])
        println(io)
    end
    nothing
end

mutable struct CanvasMatrix2x2
    A::Canvas
    B::Canvas
    C::Canvas
    D::Canvas

    nupper::Int
    nlower::Int
    nleft::Int
    nright::Int

    function CanvasMatrix2x2(
        nupper::Integer,
        nlower::Integer,
        nleft::Integer,
        nright::Integer;
        height::Union{Integer, Nothing}=nothing,
        width::Union{Integer, Nothing}=nothing,
        max_height::Union{Integer, Nothing}=nothing,
        max_width::Union{Integer, Nothing}=nothing,
        out_stream::IO=stdout,
        fix_ar::Bool=false
    )
        nrows = nupper + nlower
        ncols = nleft + nright
        height, width = get_canvas_dimensions(
            nrows,
            ncols,
            height,
            width,
            max_height,
            max_width,
            out_stream,
            fix_ar
        )
        width_A = floor(Int, width * nleft / max(ncols, 1.0))
        width_B = width - width_A
        height_A = floor(Int, height * nupper / max(nrows, 1.0))
        height_C = height - height_A
        A = Canvas(nupper, nleft, height=height_A, width=width_A)
        B = Canvas(nupper, nright, height=height_A, width=width_B)
        C = Canvas(nlower, nleft, height=height_C, width=width_A)
        D = Canvas(nlower, nright, height=height_C, width=width_B)
        new(A, B, C, D, nupper, nlower, nleft, nright)
    end
end

function point!(c::CanvasMatrix2x2, x::Integer, y::Integer)
    if x <= c.nupper
        if y <= c.nleft
            point!(c.A, x, y)
        else
            point!(c.B, x, y - c.nleft + 1)
        end
    else
        if y <= c.nleft
            point!(c.C, x - c.nupper + 1, y)
        else
            point!(c.D, x - c.nupper + 1, y - c.nleft + 1)
        end
    end
end

Base.show(io::IO, c::CanvasMatrix2x2) = _show(io, c)

function _show(io::IO, c::CanvasMatrix2x2)
    m1, m2, n1, n2 =
        size(c.A.grid, 1), size(c.C.grid, 1), size(c.A.grid, 2), size(c.B.grid, 2)
    if size(c.C.grid, 1) != size(c.D.grid, 1) || size(c.A.grid, 1) != size(c.B.grid, 1)
        @log level = 1000 "Canvas dimensions do not agree!"
    end
    for i in 1:m1
        printrow(io, c.A.grid[i, :])
        print(io, BORDER_SOLID.l)
        printrow(io, c.B.grid[i, :])
        println(io)
    end
    print(io, BORDER_SOLID.t^n1, Char(BLANK_BRAILLE), BORDER_SOLID.t^n2)
    println(io)
    for i in 1:m2
        printrow(io, c.C.grid[i, :])
        print(io, BORDER_SOLID.l)
        printrow(io, c.D.grid[i, :])
        println(io)
    end
    nothing
end
