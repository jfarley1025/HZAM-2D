using Documenter
import HZAM
makedocs(
    sitename="HZAM",
    pages=[
        "Home" => "index.md"
        "Library" => [
            "hzam.md",
            "population.md",
            "mating.md",
            "data_analysis.md",
            "plot_data.md"
        ]
        "Manual" => [
            "manual_one_simulation.md",
            "manual_simulation_set.md"
        ]
    ]
)
