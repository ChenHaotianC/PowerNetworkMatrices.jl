"""
structure and functions for power systems analysis based on Ward's decomposition
"""

"""
Data related to the considered Ward decomposition is here stored.
Given the set of buses under study, the original Ybus can be divided into 
sub matrices belongin to internal (I), boundary buses (B) and external (E) buses.
       |Y_ii  Y_ib   0  |
       |                |
 Ybus =|Y_bi  Y_bb  Y_be|
       |                |
       | 0    Y_eb  Y_ee|

# Arguments
- `y_equ_boundary::Matrix{ComplexF64}`:
        the product of Y_bb - Y_be * inv(Y_ee) * Y_eb
- `ybus_study::Matrix{ComplexF64}`:
        the whole Ybus matrix of the equavilant system
"""

struct ward_decompose{T}
    y_equ_boundary::Matrix{ComplexF64}
    ybus_study::Matrix{ComplexF64}
    location_map::Dict{Int64, Vector{Int64}}
end

function ward_decompose(
    internal_buses::Vector{Int64},
    boundary_buses::Vector{Int64},
    external_buses::Vector{Int64},
    ybus::SparseMatrixCSC{ComplexF64, Int64}
)
    # get the mapping of the index of buses
    all_buses_index = sort(vcat(internal_buses, boundary_buses, external_buses))
    bus_index_mapping = Dict(bus => index for (index, bus) in enumerate(all_buses_index))     #create a mapping for the bus index
    bus_index_inverse_mapping = Dict(index => bus for (index, bus) in enumerate(all_buses_index)) #create an inverse mapping for the bus index
    
    mapped_internal_buses = [bus_index_mapping[bus] for bus in internal_buses]
    mapped_boundary_buses = [bus_index_mapping[bus] for bus in boundary_buses]
    mapped_external_buses = [bus_index_mapping[bus] for bus in external_buses]
    
    # divide the ybus into y_ii, y_bb, y_ee,....
    y_ee = ybus[mapped_external_buses, mapped_external_buses];
    y_eb = ybus[mapped_external_buses, mapped_boundary_buses];
    y_be = ybus[mapped_boundary_buses, mapped_external_buses];
    y_bb = ybus[mapped_boundary_buses, mapped_boundary_buses];
    y_ii = ybus[mapped_internal_buses, mapped_internal_buses];
    y_ib = ybus[mapped_internal_buses, mapped_boundary_buses];
    y_bi = ybus[mapped_boundary_buses, mapped_internal_buses];
    
    # compute the equavilant boundary ybus
    y_equ_boundary = y_bb - y_be * KLU.solve!(klu(y_ee), Matrix(y_eb));
    
    # construct the ybus of study system
    study_buses = sort(union(mapped_internal_buses, mapped_boundary_buses))
    num_study_buses = length(study_buses)
    ybus_study = spzeros(ComplexF64, num_study_buses, num_study_buses) # study Ybus matrix
    # Create a mapping from bus number to index in the new matrix
    bus_to_index = Dict(bus => index for (index, bus) in enumerate(study_buses))
    
    # Create indices for boundary and internal buses in ybus_study
    boundary_indices = map(bus -> bus_to_index[bus], mapped_boundary_buses)
    internal_indices = map(bus -> bus_to_index[bus], mapped_internal_buses)
    
    # construct ybus_study with y_equ_boundary, ybus_internal, y_ib, and y_bi
    ybus_study[boundary_indices, boundary_indices] .= y_equ_boundary
    ybus_study[internal_indices, internal_indices] .= y_ii
    ybus_study[internal_indices, boundary_indices] .= y_ib
    ybus_study[boundary_indices, internal_indices] .= y_bi
    
    # relocate the external buses using minimum electrical distance
    
    num_buses = size(ybus, 1)
    lines = [(i, j, 1/ybus[i, j]) for i in 1:num_buses for j in (i + 1):num_buses if ybus[i, j] != 0.0]

    
    num_lines = length(lines)
    g = SimpleWeightedGraph(num_buses)
    for i in 1:num_lines
        line = lines[i, :]
        weight = abs(line[1][3])  # Calculate impedance magnitude
        SimpleWeightedGraphs.add_edge!(g, line[1][1], line[1][2], weight)
    end
    
    relocate_buses = copy(mapped_external_buses)
    for i in eachindex(mapped_external_buses)
        distances = Graphs.dijkstra_shortest_paths(g, mapped_external_buses[i]).dists
        min_distance, min_index = findmin(map(bus -> distances[bus], study_buses))
        relocate_buses[i] = study_buses[min_index]
    end
    
    location_map = Dict{Int64, Vector{Int64}}()  
    for (old_loc, new_loc) in zip(external_buses, relocate_buses)
        remapped_loc = bus_index_inverse_mapping[new_loc]
        push!(get!(location_map, remapped_loc, Vector{Int64}()), old_loc)
    end

    return (
        y_equ_boundary,
        ybus_study,
        location_map
        )

end