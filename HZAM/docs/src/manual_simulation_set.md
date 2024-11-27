# Running a simulation set

To run a set of simulations use the [`HZAM.run_HZAM_set`](@ref) method.

This will run the simulation for 64 different combinations of hybrid fitness and assortative mating strength. The values for hybrid fitness (`w_hyb`) are `[1, 0.95, 0.9, 0.8, 0.75, 0.6, 0.5, 0.3]` and the values for assortative mating strength (`S_AM`) are `[1, 3, 10, 30, 100, 300, 1000, Inf]`.

The outcome of each simulation is stored in the results folder which by default is `HZAM_results_GitIgnore/simulation_outcomes`. The results folder can be changed using the [`HZAM.set_results_folder`](@ref) method.

[`HZAM.run_HZAM_set`](@ref) returns an outcome array with all the data from each simulation. This array is not automatically saved to file (though it contains the same information that's stored in the simulation set output folder). To save it to a file use the JLD2 package like so:

```julia
using JLD2
import HZAM

outcome_array = HZAM.run_HZAM_set("my_set", 1.1, 0)

@save "my_file.JLD2" outcome_array
```

## Working with the output data

The [`HZAM.run_HZAM_set`](@ref) automatically saves the outcome of each simulation to a file. These files can either be read individually or the entire folder can be loaded into an outcome array using [`HZAM.load_from_folder`](@ref).

Data can also be converted to CSV format, but it first needs to be extracted from the OutputData structure.

Here's an example of how to do this:

```julia
import HZAM

outcome_array = HZAM.load_from_folder("HZAM_results_GitIgnore/simulation_outcomes/myset")

bimodality_array = [o.spatial_data.bimodality for o in outcome_array]
w_hyb_array = [o.sim_params.w_hyb for o in outcome_array]
S_AM_array = [o.sim_params.S_AM for o in outcome_array]

HZAM.onvert_to_CSV(
    bimodality_array,
    w_hyb_array,
    S_AM_array,
    "bimodality"
)
```



