using Revise
using Test
import Logging
using PowerNetworkMatrices
using PowerSystems
using InfrastructureSystems
using PowerSystemCaseBuilder
using TimeSeries
using LinearAlgebra
using Plots
using BenchmarkTools

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

# get systems
sys2 = System("ACTIVSg2000.m")

@benchmark begin
    # get the PTDF matrix
    ptdf = PTDF(sys)
    a = IncidenceMatrix(sys)

    # get covariance matrix --> it is square
    ptdf_1 = ptdf.data[:, setdiff(1:end, a.ref_bus_positions)]

    # get the norm for each row and anlge between rows
    norm_ptdf = [norm(ptdf_1[i, :]) for i in axes(ptdf_1, 1)]
    angles = Array{Array{Float64}}(undef, size(ptdf_1, 1) - 1)
    for i in 1:(size(ptdf_1, 1) - 1)
        angles[i] =
            acosd.(
                clamp.(
                    (ptdf_1[(i + 1):end, :] * ptdf_1[i, :]) ./
                    (norm_ptdf[(i + 1):end] .* norm_ptdf[i]),
                    -1, 1)
            )
    end
end

"""
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 10.466 s (8.13% GC) to evaluate,
 with a memory estimate of 76.94 GiB, over 116343 allocations.
"""

# get all angles
all_angles = sort(reduce(vcat, angles))
for i in 1:length(angles)
    if sum(angles[i] .< 10) > 3 || sum(angles[i] .> 170) > 3
        @show i
    end
end

# sort and plot
plot(all_angles)
histogram(all_angles)
