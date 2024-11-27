
using Distributed
interrupt()

@time begin
	using Pkg
	Pkg.activate("HZAM 2.0/HZAM")
	Pkg.instantiate()
	Pkg.precompile()
end

rmprocs(20)
addprocs(4)

@everywhere begin
	using Pkg
	Pkg.activate("HZAM 2.0/HZAM")
	Pkg.instantiate()
	Pkg.precompile()
end

@everywhere begin
	include("HZAM/src/HZAM.jl")

	import .HZAM
	using JLD2 # needed for saving / loading data in Julia format
	HZAM.set_results_folder(string("HZAM-J_2D_results/simulation_outcomes/Run2_2024July22"))
end

# run main simulation sets
HZAM.run_HZAM_sets_complete("Run2_20240722")
