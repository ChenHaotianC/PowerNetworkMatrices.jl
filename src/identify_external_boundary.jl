function identify_external_boundary(
    all_lines::Vector{ACBranch},
    internal_buses::Vector{Int64}
)
    from_buses = [line.arc.from.number for line in all_lines]
    to_buses = [line.arc.to.number for line in all_lines]

    all_buses = unique(union(from_buses,to_buses));

    # Identify boundary buses: those connected to internal buses
    is_from_internal = [bus in internal_buses for bus in from_buses]
    is_to_internal = [bus in internal_buses for bus in to_buses]
    boundary_from_buses = [to_buses[i] for i in 1:length(from_buses) if is_from_internal[i] && !is_to_internal[i]]
    boundary_to_buses = [from_buses[i] for i in 1:length(to_buses) if !is_from_internal[i] && is_to_internal[i]]
    boundary_buses = unique(vcat(boundary_from_buses, boundary_to_buses))

    # Identify external buses
    is_external = .!(in.(all_buses, Ref(internal_buses)).| in.(all_buses, Ref(boundary_buses)))
    external_buses = all_buses[is_external]

    return (boundary_buses,external_buses)
end