function legendre_(order::Int, x::Number)
    @assert(order >= 0, "order cannot be <0")
    @assert(order < 11, "order cannot be <0")
    order == 0 && return 1
    order == 1 && return x
    order == 2 && return 0.5 * (3 * x^2 - 1)
    order == 3 && return 0.5 * (5 * x^3 - 3 * x)
    order == 4 && return (1.0 / 8) * (35 * x^4 - 30 * x^2 + 3)
    order == 5 && return (1.0 / 8) * (63 * x^5 - 70 * x^3 + 15 * x)
    order == 6 && return (1.0 / 16) * (
        231 * x^6 - 315 * x^4 + 105 * x^2 - 5
    )
    order == 7 && return (1.0 / 16) * (
        429 * x^7 - 693 * x^5 + 315 * x^3 - 35 * x
    )
    order == 8 && return (1.0 / 128) * (
        6435 * x^8
        -
        12012 * x^6
        +
        6930 * x^4
        -
        1260 * x^2
        +
        35
    )
    order == 9 && return (1.0 / 128) * (
        12155 * x^9
        -
        25740 * x^7
        +
        18018 * x^5
        -
        4620 * x^3
        +
        315 * x
    )
    order == 10 && return (1.0 / 256) * (
        46189 * x^10
        -
        109395 * x^8
        +
        90090 * x^6
        -
        30030 * x^4
        +
        3465 * x^2
        -
        63
    )
end



function legendre(order::Int, x::Number)
    return sqrt(2 * order + 1) * legendre_(order, x)
end


function local_to_global_coords(nc::Int, cell::Int, local_coord::Number)
    @assert(local_coord >= -1, "local_coord >= -1")
    @assert(local_coord <= 1, "local_coord <= 1")
    return (local_coord / 2 + 0.5) * ((1 + cell) / nc)
end

function global_to_local_coords(x_global, c, nc, cell_size)
    return (((x_global - c * cell_size) / cell_size) - 1) * 2 + 1
end
