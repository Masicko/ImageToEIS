module ImageToEIS

#using ElectricalEquivalentCircuits
using ImageIO
using FileIO
using IterativeSolvers
using PyPlot
using Dates
using Printf
using Colors
using DataFrames
using CSV
using Interpolations
using Optim
using SparseArrays

using Base

include("physical_properties.jl")
export matrix_to_file

# general part

include("material_matrix_to_lin_sys.jl")
include("material_matrix_to_impedance.jl")
export material_matrix_to_impedance

# user requested part

include("equivalent_circuit_support.jl")
include("Z_view_export.jl")
include("image_to_EIS_interface.jl")
export image_to_EIS

# maybe independent package for proper geometry generation

include("generate_matrix.jl")
export generate_matrix
export three_column_domain_template
export three_column_domain_matrix
export chess_matrix

# macro stuff consisting of par study, evaluate results, ploting, cluster operations

include("macro_stuff.jl")
export par_study

# specific material dependences

include("temperature_dependences.jl")
export TI


end # module
