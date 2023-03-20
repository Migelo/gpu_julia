function get_polynomal_order(k::Int)::Int
    k + 1
end

function get_number_of_internal_qpoints(number_of_dimensions::Int, number_of_qpoints_per_dimension)::Int
    number_of_qpoints_per_dimension^number_of_dimensions
end

function get_number_of_outer_qpoints(number_of_dimensions::Int, number_of_qpoints_per_dimension)::Int
    number_of_qpoints_per_dimension^(number_of_dimensions - 1)
end

function get_number_of_runge_kutta_stages(dg_order::Int)::Int
    if dg_order < 4
        return dg_order
    else
        return 5
    end
end

function get_number_of_cells(number_of_dimensions::Int, cells_per_dimension::Int)::Int
    return cells_per_dimension^number_of_dimensions
end

function get_total_number_of_weights(number_of_cells::Int, number_of_fields::Int, number_of_orders::Int)::Int
    number_of_cells * number_of_fields * number_of_orders
end


DEGREE_K = 1
ND = 1
CELLS_PER_DIMENSION = 32
NF = (1 + ND)

NB = get_polynomal_order(DEGREE_K)
QPOINTS_PER_DIMENSION = NB
NUM_INNER_QPOINTS = get_number_of_internal_qpoints(ND, QPOINTS_PER_DIMENSION)
NUM_OUTER_QPOINTS = get_number_of_outer_qpoints(ND, QPOINTS_PER_DIMENSION)
RK_STAGES = get_number_of_runge_kutta_stages(DEGREE_K)
NC = get_number_of_cells(ND, CELLS_PER_DIMENSION)
NW = get_total_number_of_weights(NC, NF, NB)

