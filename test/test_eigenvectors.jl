using Revise
using Test
import Logging
using PowerNetworkMatrices
using PowerSystems
using InfrastructureSystems
using PowerSystemCaseBuilder
using TimeSeries
using LinearAlgebra

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

# get systems
sys = System("ACTIVSg2000.m")

# get the PTDF matrix
ptdf = PTDF(sys)
a = IncidenceMatrix(sys)

# the matrix rankk is equal to the min(n_rows, n_columns - n_slack_buses)
rank(ptdf.data[:, setdiff(1:end, a.ref_bus_positions)]) == size(ptdf.data[:, setdiff(1:end, a.ref_bus_positions)], 2)
# get covariance matrix --> it is square
ptdf_1 = ptdf.data[:, setdiff(1:end, a.ref_bus_positions)]
co_ptdf_1 = ptdf_1*ptdf_1'

# get eigs and eigenvects of covarianc matrix 
# --> they indicate the direcation where data changes the most
eigenvals = eigvals(co_ptdf_1)
eigenvects = eigvecs(co_ptdf_1)

# get the basis of the ptdf matrix --> linearly independent vectors

# evaluate angle between them

# evaluate angle between all rows