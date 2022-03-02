module ImageToEIS

#using ElectricalEquivalentCircuits
using ImageIO
using FileIO
using IterativeSolvers
using PyPlot
using Printf

using Base

include("auxilary_stuff.jl")
include("Z_view_export.jl")
include("physical_properties.jl")
include("image_to_equation_matrix.jl")
include("impedance_evaluation.jl")
export image_to_EIS

include("macro_stuff.jl")




end # module
