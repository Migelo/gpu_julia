function legendre(order::Int, x::Number)::Float64
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



function legendre_normalized(order::Int, x::Number)
    return sqrt(2 * order + 1) * legendre(order, x)
end

function legendre_diff(order::Int, x::Float64)::Float64
    if order == 0
        return 0.0
    elseif order == 1
        return 1.0
    elseif order == 2
        return 3.0 * x
    elseif order == 3
        return 0.5 * (15 * x^2 - 3)
    elseif order == 4
        return (1.0 / 8) * (35 * 4 * x^3 - 60 * x)
    elseif order == 5
        return (1.0 / 8) * (63 * 5 * x^4 - 70 * 3 * x^2 + 15)
    elseif order == 6
        return (1.0 / 16) * (231 * 6 * x^5 - 315 * 4 * x^3 + 105 * 2 * x)
    elseif order == 7
        return (1.0 / 16) * (429 * 7 * x^6 - 693 * 5 * x^4 + 315 * 3 * x^2 - 35)
    elseif order == 8
        return (1.0 / 128) * (6435 * 8 * x^7 - 12012 * 6 * x^5 + 6930 * 4 * x^3 - 1260 * 2 * x)
    elseif order == 9
        return (1.0 / 128) * (12155 * 9 * x^8 - 25740 * 7 * x^6 + 18018 * 5 * x^4 - 4620 * 3 * x^2 + 315)
    elseif order == 10
        return (1.0 / 256) * (46189 * 10 * x^9 - 109395 * 8 * x^7 + 90090 * 6 * x^5 - 30030 * 4 * x^3 + 3465 * 2 * x)
    else
        error("undefined")
    end
end



function local_to_global_coords(nc::Int, cell::Int, local_coord::Number)::Float64
    @assert(local_coord >= -1, "local_coord >= -1")
    @assert(local_coord <= 1, "local_coord <= 1")
    return (local_coord / 2 + 0.5) * ((1 + cell) / nc)
end

function global_to_local_coords(x_global, c, cell_size)::Float64
    return (((x_global - c * cell_size) / cell_size) - 1) * 2 + 1
end
