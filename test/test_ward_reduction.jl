@testset "Test ward reduction" begin
    # test on 5, 14, 30 and 73 bus system
    for name in ["c_sys5","c_sys14","pti_case73_sys","pti_case30_sys"]
        if startswith(name, "c")
            sys = PSB.build_system(PSB.PSITestSystems, name)
        elseif startswith(name, "pti")
            sys = PSB.build_system(PSB.PSSEParsingTestSystems, name)
        end
        # get the lines data
        all_lines = get_components(ACBranch, sys);
        all_lines = collect(all_lines)

        # get the mapping index
        all_buses_index = PSY.get_bus_numbers(sys)
        bus_index_mapping = Dict(bus => index for (index, bus) in enumerate(all_buses_index))  #create a mapping for the bus index
        inverse_bus_index_mapping = Dict(index => bus for (index, bus) in enumerate(all_buses_index))
        # get the ybus matrix from system
        ybus = PNM.Ybus(sys).data
        # get the PTDF from the system
        ptdf = PNM.PTDF(sys).data
        # get the reference bus
        all_buses = get_components(ACBus, sys)
        ref_bus = [bus for bus in all_buses if bus.bustype == ACBusTypes.REF][1].number
        mapped_ref_bus = bus_index_mapping[ref_bus]; 
        # compute the ptdf from ybus
        ptdf_compute = ptdf_by_ybus(Matrix(ybus),mapped_ref_bus)
        # test if the PTDF computed by Ybus is the same as the ptdf from system
        number_buses = size(ptdf,1);
        number_lines = size(ptdf,2);
        for i in 1:number_buses
            for j in 1:number_lines
                # @test isapprox(ptdf[i,j], transpose(ptdf_compute)[i, j], atol = 1e-4)
            end
        end
        
        # drop the resistance and make the ybus matrix singular
        ybus_original = map(x -> x != 0.0 ? -1/imag(1/x)*im : x, ybus)
        row_sums = sum(ybus_original, dims=2) - diag(ybus_original)
        ybus_original[diagind(ybus_original)] = -row_sums

        # select the internal buses
        pick_zone = '1';
        internal_buses = [index for index in all_buses_index if string(index)[1] == pick_zone];
        boundary_buses,external_buses = identify_external_boundary(all_lines,internal_buses);

        # get the mapping internal boundary and external buses
        mapped_internal_buses = [bus_index_mapping[bus] for bus in internal_buses]
        mapped_boundary_buses = [bus_index_mapping[bus] for bus in boundary_buses]
        mapped_external_buses = [bus_index_mapping[bus] for bus in external_buses]

        ptdf_compute_original = ptdf_by_ybus(Matrix(ybus_original),mapped_ref_bus)
        ward_data = ward_decompose(internal_buses,boundary_buses,external_buses,ybus_original);
        # ward_data = @btime ward_decompose($internal_buses, $boundary_buses, $external_buses, $ybus_original) # check the time of ward decompose
        ybus_study = ward_data[2] # the ybus of study system

        # set the new reference bus
        num_less = count(x -> x < mapped_ref_bus, mapped_external_buses)
        ref_bus_study = mapped_ref_bus - num_less
        ptdf_study = ptdf_by_ybus(Matrix(ybus_study),ref_bus_study)

        # find the corresponding index of buses and branches in the original system and study system
        # get the index of internal branches in original system
        from_to_buses_original = [(inverse_bus_index_mapping[i], inverse_bus_index_mapping[j]) for i in 1:number_buses for j in (i+1):number_buses if ybus_original[i, j] != 0.0]
        is_from_internal = in.(first.(from_to_buses_original), Ref(internal_buses))
        is_to_internal = in.(last.(from_to_buses_original), Ref(internal_buses))
        is_from_external = in.(first.(from_to_buses_original), Ref(external_buses))
        is_to_external = in.(last.(from_to_buses_original), Ref(external_buses))
        condition1 = is_from_internal .& .!is_to_external
        condition2 = is_to_internal .& .!is_from_external
        indices_connecting_internal = findall(condition1 .| condition2)

        # get the index of internal branches in study system
        combined_buses = union(internal_buses, boundary_buses)
        sorted_buses = sort(combined_buses)
        bus_index_mapping_study = Dict(i => sorted_buses[i] for i in eachindex(sorted_buses))

        from_to_buses_study = [(bus_index_mapping_study[i], bus_index_mapping_study[j]) for i in axes(ybus_study, 1) for j in (i+1):size(ybus_study, 2) if ybus_study[i, j] != 0.0]
        from_bus_mask = [from in internal_buses for (from, to) in from_to_buses_study]
        to_bus_mask = [to in internal_buses for (from, to) in from_to_buses_study]
        combined_mask = from_bus_mask .| to_bus_mask
        indices_connecting_internal_study = findall(combined_mask)

        i_b_buses = sort(union(mapped_internal_buses,mapped_boundary_buses))
        sub_ptdf_original = ptdf_compute_original[indices_connecting_internal, i_b_buses]
        sub_ptdf_study = ptdf_study[indices_connecting_internal_study, 1:length(i_b_buses)]
        @test all(isapprox.(sub_ptdf_original, sub_ptdf_study, atol=1e-4))  # compare if the corresponding elements in PTDFs are the same
    end
end