"Types and functions used for analyzing data produced by the simulation."
module DataAnalysis

using LsqFit: curve_fit

mean(itr) = sum(itr) / length(itr)

include("used_data_analysis.jl")

"unused_data_analysis.jl contains some extra functions that are not currently in use."

# include("unused_data_analysis.jl")
end