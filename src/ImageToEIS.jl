module ImageToEIS

#using ElectricalEquivalentCircuits
using ImageIO
using FileIO
using IterativeSolvers
using PyPlot
using Dates
using Printf

using Base

include("equivalent_circuit_support.jl")
include("Z_view_export.jl")
include("physical_properties.jl")
include("image_to_equation_matrix.jl")
include("matrix_to_impedance.jl")
include("image_to_EIS_interface.jl")
export image_to_EIS

include("macro_stuff.jl")




end # module
