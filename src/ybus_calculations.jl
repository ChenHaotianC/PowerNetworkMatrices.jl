"""
Nodal admittance matrix (Ybus) is an N x N matrix describing a power system with N buses. It represents the nodal admittance of the buses in a power system.

The Ybus Struct is indexed using the Bus Numbers, no need for them to be sequential
"""
struct Ybus{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF64}
    data::SparseArrays.SparseMatrixCSC{ComplexF64, Int}
    axes::Ax
    lookup::L
end

function _ybus!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    br::PSY.ACBranch,
    num_bus::Dict{Int, Int},
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]
    Y_l = (1 / (PSY.get_r(br) + PSY.get_x(br) * 1im))
    Y11 = Y_l + (1im * PSY.get_b(br).from)
    if !isfinite(Y11) || !isfinite(Y_l)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end
    ybus[bus_from_no, bus_from_no] += Y11
    Y12 = -Y_l
    ybus[bus_from_no, bus_to_no] += Y12
    #Y21 = Y12
    ybus[bus_to_no, bus_from_no] += Y12
    Y22 = Y_l + (1im * PSY.get_b(br).to)
    ybus[bus_to_no, bus_to_no] += Y22
    return
end

function _ybus!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    br::PSY.DynamicBranch,
    num_bus::Dict{Int, Int},
)
    _ybus!(ybus, br.branch, num_bus)
    return
end

function _ybus!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    br::PSY.Transformer2W,
    num_bus::Dict{Int, Int},
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    Y11 = Y_t
    b = PSY.get_primary_shunt(br)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(b)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    ybus[bus_from_no, bus_from_no] += Y11 - (1im * b)
    ybus[bus_from_no, bus_to_no] += -Y_t
    ybus[bus_to_no, bus_from_no] += -Y_t
    ybus[bus_to_no, bus_to_no] += Y_t
    return
end

function _ybus!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    br::PSY.TapTransformer,
    num_bus::Dict{Int, Int},
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]

    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    c = 1 / PSY.get_tap(br)
    b = PSY.get_primary_shunt(br)

    Y11 = (Y_t * c^2)
    ybus[bus_from_no, bus_from_no] += Y11 - (1im * b)
    Y12 = (-Y_t * c)
    if !isfinite(Y11) || !isfinite(Y12) || !isfinite(b)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end
    ybus[bus_from_no, bus_to_no] += Y12
    #Y21 = Y12
    ybus[bus_to_no, bus_from_no] += Y12
    Y22 = Y_t
    ybus[bus_to_no, bus_to_no] += Y22
    return
end

function _ybus!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    br::PSY.PhaseShiftingTransformer,
    num_bus::Dict{Int, Int},
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    tap = (PSY.get_tap(br) * exp(PSY.get_α(br) * 1im))
    c_tap = (PSY.get_tap(br) * exp(-1 * PSY.get_α(br) * 1im))
    b = PSY.get_primary_shunt(br)
    Y11 = (Y_t / abs(tap)^2)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(b * c_tap)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end
    ybus[bus_from_no, bus_from_no] += Y11 - (1im * b)
    Y12 = (-Y_t / c_tap)
    ybus[bus_from_no, bus_to_no] += Y12
    Y21 = (-Y_t / tap)
    ybus[bus_to_no, bus_from_no] += Y21
    Y22 = Y_t
    ybus[bus_to_no, bus_to_no] += Y22
    return
end

function _ybus!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    fa::PSY.FixedAdmittance,
    num_bus::Dict{Int, Int},
)
    bus = PSY.get_bus(fa)
    bus_no = num_bus[PSY.get_number(bus)]
    if !isfinite(fa.Y)
        error(
            "Data in $(PSY.get_name(fa)) is incorrect. Y = $(fa.Y)",
        )
    end
    ybus[bus_no, bus_no] += fa.Y
    return
end

function _buildybus(
    branches,
    buses::Vector{PSY.Bus},
    fixed_admittances::Vector{PSY.FixedAdmittance},
)
    buscount = length(buses)
    num_bus = Dict{Int, Int}()

    for (ix, b) in enumerate(buses)
        num_bus[PSY.get_number(b)] = ix
    end
    ybus = SparseArrays.spzeros(ComplexF64, buscount, buscount)
    for b in branches
        if PSY.get_name(b) == "init"
            throw(DataFormatError("The data in Branch is invalid"))
        end
        PSY.get_available(b) && _ybus!(ybus, b, num_bus)
    end
    for fa in fixed_admittances
        PSY.get_available(fa) && _ybus!(ybus, fa, num_bus)
    end
    return SparseArrays.dropzeros!(ybus)
end

"""
Builds a Ybus from a collection of buses and branches. The return is a Ybus Array indexed with the bus numbers and the branch names.

# Arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Depth First Search (DFS)
"""
function Ybus(
    branches::Vector,
    buses::Vector{PSY.Bus},
    fixed_admittances::Vector{PSY.FixedAdmittance} = Vector{PSY.FixedAdmittance}();
    check_connectivity::Bool = true,
)
    bus_ax = PSY.get_number.(buses)
    axes = (bus_ax, bus_ax)
    bus_lookup = make_ax_ref(bus_ax)
    look_up = (bus_lookup, bus_lookup)
    ybus = _buildybus(branches, buses, fixed_admittances)
    if check_connectivity && length(buses) > 1
        islands = find_subnetworks(ybus, bus_ax)
        length(islands) > 1 && throw(IS.DataFormatError("Network not connected"))
    end
    return Ybus(ybus, axes, look_up)
end

"""
Builds a Ybus from the system. The return is a Ybus Array indexed with the bus numbers and the branch names.

# Arguments
- `check_connectivity::Bool`: Checks connectivity of the network
"""
function Ybus(sys::PSY.System; kwargs...)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    fixed_admittances = collect(PSY.get_components(PSY.FixedAdmittance, sys))
    return Ybus(
        branches,
        buses,
        fixed_admittances;
        kwargs...,
    )
end

"""
Finds the equivalent impedance between nodes in a network.

# Arguments
- `ybus::Ybus.
"""
function find_equivalent_impedance(ybus::Ybus, from::Int, to::Int)



function _find_equivalent_impedance(M::SparseArrays.SparseMatrixCSC{ComplexF64, Int}, bus_numbers::Vector{Int})
    rows = SparseArrays.rowvals(M)
    touched = Set{Int}()
    subnetworks = Dict{Int, Set{Int}}()
    for (ix, bus_number) in enumerate(bus_numbers)
        neighbors = SparseArrays.nzrange(M, ix)
        if length(neighbors) < 1
            @warn "Bus $bus_number is islanded"
            subnetworks[bus_number] = Set{Int}(bus_number)
            continue
        end
        for j in neighbors
            row_ix = rows[j]
            if bus_number ∉ touched
                push!(touched, bus_number)
                subnetworks[bus_number] = Set{Int}(bus_number)
                _dfs(row_ix, M, bus_numbers, subnetworks[bus_number], touched)
            end
        end
    end
    return subnetworks
end
