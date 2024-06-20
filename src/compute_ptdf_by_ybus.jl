function ptdf_by_ybus(
    ybus::Matrix{ComplexF64},
    ref_bus::Int64,
)
    num_buses = size(ybus, 1)

    # Identify branches from non-zero off-diagonal elements of Ybus
    lines = [(i, j) for i in 1:num_buses for j in (i+1):num_buses if ybus[i, j] != 0.0]
    num_lines = length(lines)

    # Incidence Matrix (A)
    Incidence_Matrix = zeros(Int, num_buses, num_lines)
    for (line_idx, (from_bus, to_bus)) in enumerate(lines)
        Incidence_Matrix[from_bus, line_idx] = 1
        Incidence_Matrix[to_bus, line_idx] = -1
    end

    Incidence_Matrix = transpose(Incidence_Matrix)

    # Branch Admittance Matrix (BA_matrix)
    BA_matrix = zeros(ComplexF64, num_lines, num_buses)
    for (line_idx, (from_bus, to_bus)) in enumerate(lines)
        BA_matrix[line_idx, from_bus] = 1/imag(1/ybus[from_bus, to_bus])
        BA_matrix[line_idx, to_bus] = -1/imag(1/ybus[from_bus, to_bus])
    end
    BA_matrix = real(BA_matrix)

    # compute the Bbus matrix
    Bbus = transpose(Incidence_Matrix) * BA_matrix

    # compute the PTDF
    H = zeros(num_lines, num_buses)
    Bf_noref = BA_matrix[:, 1:end .!= ref_bus]  # Remove the column corresponding to the reference bus
    Bbus_noref = Bbus[1:end .!= ref_bus, 1:end .!= ref_bus]  # Remove both row and column for the reference bus
    H[:, 1:end .!= ref_bus] = Bf_noref / Bbus_noref  # Using matrix division for the inversion and multiplication
    H[:, ref_bus] = zeros(size(H, 1))

    return H
end