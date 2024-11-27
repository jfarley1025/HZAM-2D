include("HZAM/src/HZAM.jl")

import .HZAM
using JLD2 # needed for saving / loading data in Julia format

using Plots
using BenchmarkTools

outcome =
    HZAM.run_one_HZAM_sim(1, 1000, 1.1; # these values are 
        # hybrid fitness; AM strength; ecol. diff; intrinsic growth rate 
        K_total=30000, max_generations=1000,
        do_plot=false, plot_int=10,
        total_loci=1,
        female_mating_trait_loci=[1],
        male_mating_trait_loci=[1],
        hybrid_survival_loci=[1],
        per_reject_cost=0,
        sigma_disp=0.03,
        sigma_comp=0.01,
        track_population_data=true,
        exit_early = false
    )

    @save "sample_outcome.jld2" outcome