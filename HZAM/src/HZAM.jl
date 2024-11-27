#=import Pkg
Pkg.activate("./HZAM")
=#

module HZAM
include("data_analysis.jl") # functions for analyzing data (calculating width/length of hybrid zone, gene flow, bimodality, etc.)
include("population.jl") # data types and functions for initializing population with genotypes, locations, growth rates, etc.
include("mating.jl") # functions for finding a mate and determining match strength
include("plot_data.jl")

import .DataAnalysis
import .Mating
import .Population.NUM_ZONES
import .PlotData
import .Population
using Distributed
using JLD2

include("simulation.jl") # main function for running a simulation
include("summarize_data.jl")

end # module HZAM
