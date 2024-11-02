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

@testset "Test ward reduction" begin
    sys = PSB.build_system(PSB.PSSEParsingTestSystems, "pti_case73_sys");
    all_areas = PSY.get_components(PSY.Area, sys);
    for area in collect(all_areas)

        #identify the internal buses which are located in the selected area
        internal_buses = sort([parse(Int64, get_name(b)) for b in PSY.get_components(PSY.Bus, sys) if get_area(b) == area])

        # use the ward decompose to get the ybus matrix of study system
        ybus_study = ward_decompose(sys,internal_buses)                 # the decomposed ybus
        ybus_study_matrix = ybus_study.data;

        # for testing, compare the elements of PTDF computed by ybus_original and ybus_study

        # compute the PTDF based on the  ybus
        ybus = PNM.Ybus(sys)
        ybus_matrix = ybus.data;
        # drop the resistance and make the ybus matrix singular
        ybus_matrix = map(x -> x != 0.0 ? -1/imag(1/x)*im : x, ybus_matrix)
        row_sums = sum(ybus_matrix, dims=2) - diag(ybus_matrix)
        ybus_matrix[diagind(ybus_matrix)] = -row_sums
        
        buses = PSY.get_bus_numbers(sys)
        bus_lookup = PNM.make_ax_ref(buses) 
        ref_bus = internal_buses[1]                           # select the first bus as the reference bus

        ptdf_original = ptdf_by_ybus(Matrix(ybus_matrix),bus_lookup[ref_bus])                            #the PTDF of original system
        ptdf_study = ptdf_by_ybus(Matrix(ybus_study_matrix),ybus_study.lookup[1][ref_bus])        #the PTDF of study system

        #compare the ptdf of corresponding buses and line
        ybus_original_lookup = ybus.lookup[1]
        ybus_study_lookup = ybus_study.lookup[1]
        original_branches = []
        for i in axes(ybus_matrix, 1)
            for j in (i+1):size(ybus_matrix, 1)
                if ybus_matrix[i, j] != 0
                    push!(original_branches, (first([k for (k, v) in ybus_original_lookup if v == i]), first([k for (k, v) in ybus_original_lookup if v == j])))
                end
            end
        end      
        study_branches = []
        for i in axes(ybus_study_matrix, 1)
            for j in (i+1):size(ybus_study_matrix, 1)
                if ybus_study_matrix[i, j] != 0
                    push!(study_branches, (first([k for (k, v) in ybus_study_lookup if v == i]), first([k for (k, v) in ybus_study_lookup if v == j])))
                end
            end
        end      
        common_branches = intersect(original_branches, study_branches)
        original_indices = [first(findall(x -> x == branch, original_branches)) for branch in [(from, to) for (from, to) in common_branches if from in internal_buses && to in internal_buses]]
        study_indices = [first(findall(x -> x == branch, study_branches)) for branch in [(from, to) for (from, to) in common_branches if from in internal_buses && to in internal_buses]]
        
        for bus_number in buses
            if haskey(ybus_original_lookup, bus_number) && haskey(ybus_study_lookup, bus_number)
               @test all(isapprox.(ptdf_original[original_indices, ybus_original_lookup[bus_number]], ptdf_study[study_indices, ybus_study_lookup[bus_number]], atol=1e-4))
            end
        end
    end
end