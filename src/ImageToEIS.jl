module ImageToEIS

#using ElectricalEquivalentCircuits
using ImageIO
using FileIO
using IterativeSolvers
using PyPlot
using Dates
using Printf
using Colors

using Base

include("equivalent_circuit_support.jl")
include("Z_view_export.jl")
include("physical_properties.jl")
export matrix_to_file

include("image_to_equation_matrix.jl")
include("matrix_to_impedance.jl")
include("image_to_EIS_interface.jl")
export image_to_EIS

include("generate_matrix.jl")
export generate_matrix
export three_column_domain_template

include("macro_stuff.jl")




end # module
