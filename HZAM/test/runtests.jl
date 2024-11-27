using Test
import HZAM.Population as Population
using HZAM.Population
import HZAM.Mating as Mating
import HZAM.DataAnalysis as DataAnalysis

test_files = [
    "test_data_analysis.jl",
    "test_find_mate.jl",
    "test_growth_rate_calculation.jl",
    "test_population.jl"
]

if length(ARGS) > 0 && !isempty(ARGS[1])
    filtered_test_files = []
    for test_file in test_files
        if test_file == ARGS[1]
            push!(filtered_test_files, test_file)
        end
    end
    test_files = filtered_test_files
end

for fname in test_files
    println(fname)
    include(fname)
end