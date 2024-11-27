using Colors, ColorSchemes
#using GLMakie
using JLD2 # needed for saving / loading data in Julia format
using LsqFit: curve_fit
using Statistics

"Compute the mean value of a vector."
mean(itr) = sum(itr) / length(itr)

"The default directory where all files are saved."
global results_folder = mkpath("HZAM-J_2D_results")

"The active directory that simulation.jl uses."
global working_dir = results_folder


"The set of hybrid fitnesses (w_hyb) values that will be run"
global w_hyb_set = [1, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]


"The set of assortative mating strengths (S_AM) values that will be run"
global S_AM_set = [1, 3, 10, 30, 100, 300, 1000, Inf]

"Evenly spaced x coordinates over the range."
global spaced_locations = collect(Float32, 0:0.01:1)

"""
	run_HZAM_sets_complete(run_name::String;
		set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
		w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
		S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
		max_generations::Integer = 2000,
		K_total::Integer = 30000,
	)

Run the main set of simulations used in the paper and save the outcomes to `results_folder`.

# Arguments

- `run_name::String`: the name of the run.
- `set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: which sets to run (the sets 
are full pleiotropy, no pleiotropy, magic cue, magic preference, etc.).
- `w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which w_hyb values to use.
- `S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which S_AM values to use.
- `max_generations::Integer`: how many generations to run the simulations for.
- `K_total::Integer`: the carrying capacity of the simulated range.
"""
function run_HZAM_sets_complete(run_name::String;
	set_numbers::Union{UnitRange{<:Integer}, Vector{<:Integer}} = 1:9,
	w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
	S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
	max_generations::Integer = 2000,
	K_total::Integer = 30000,
)
	set_names = ["full_pleiotropy", "no_pleiotropy", "separate_mmt", "separate_fmt",
		"separate_hst", "low_reject_full_pleiotropy", "high_reject_full_pleiotropy",
		"low_reject_no_pleiotropy", "high_reject_no_pleiotropy",
		"one_locus_full_pleiotropy", "one_locus_no_pleiotropy"]
	total_loci = [6, 12, 9, 9, 9, 6, 6, 9, 9, 1, 3]
	female_mating_trait_loci = [1:3, 4:6, 1:3, 4:6, 4:6, 1:3, 1:3, 4:6, 4:6, [1], [1]]
	male_mating_trait_loci = [1:3, 7:9, 4:6, 1:3, 4:6, 1:3, 1:3, 7:9, 7:9, [1], [2]]
	hybrid_survival_loci = [1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, [1], [3]]
	per_reject_cost = [0, 0, 0, 0, 0, 0.01, 0.05, 0.01, 0.05, 0, 0]

	println("$(nprocs()) threads running")
	@async @distributed for i in set_numbers
		run_HZAM_set(
			set_names[i],
			total_loci[i],
			female_mating_trait_loci[i],
			male_mating_trait_loci[i],
			hybrid_survival_loci[i],
			per_reject_cost[i],
			w_hyb_set_of_run,
			S_AM_set_of_run;
			max_generations,
			K_total,
		)
		println("--------------------")
		println(string(set_names[i], " completed successfully!"))
		println("--------------------")
	end
end

"""
	run_HZAM_sets_supplemental(run_name::String;
		w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
		S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
		max_generations::Integer = 2000,
		K_total::Integer = 30000,
	)

Run an extra set of simulations and save the outcomes to `results_folder`. Currently setup 
to run no pleiotropy and full pleitropy sets with one locus.

# Arguments

- `run_name::String`: the name of the run.
- `w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which w_hyb values to use.
- `S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}}`: which S_AM values to use.
- `max_generations::Integer`: how many generations to run the simulations for.
- `K_total::Integer`: the carrying capacity of the simulated range.
"""
function run_HZAM_sets_supplemental(run_name::String;
	w_hyb_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = w_hyb_set,
	S_AM_set_of_run::Union{UnitRange{<:Real}, Vector{<:Real}} = S_AM_set,
	max_generations::Integer = 2000,
	K_total::Integer = 30000,
)
	set_names = ["one_locus_full_pleiotropy", "one_locus_no_pleiotropy",
		"two_loci_full_pleiotropy", "two_loci_no_pleiotropy"]
	total_loci = [1, 3, 2, 6]
	female_mating_trait_loci = [[1], [2], 1:2, 3:4]
	male_mating_trait_loci = [[1], [3], 1:2, 5:6]
	hybrid_survival_loci = [[1], [1], 1:2, 1:2]
	per_reject_cost = [0, 0, 0, 0]

	println("$(nprocs()) threads running")
	@async @distributed for i in set_numbers
		run_HZAM_set(
			set_names[i],
			total_loci[i],
			female_mating_trait_loci[i],
			male_mating_trait_loci[i],
			hybrid_survival_loci[i],
			per_reject_cost[i],
			w_hyb_set_of_run,
			S_AM_set_of_run;
			max_generations,
			K_total,
		)
		println("--------------------")
		println(string(set_names[i], " completed successfully!"))
		println("--------------------")
	end
end


"""
	run_HZAM_set(
		set_name::String,
		total_loci::Int,
		female_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		male_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		hybrid_survival_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		per_reject_cost::Real,
		w_hyb_set::Union{UnitRange{<:Real}, Vector{<:Real}},
		S_AM_set::Union{UnitRange{<:Real}, Vector{<:Real}};
		<keyword arguments>
	)

Run the simulation for the given combinations of hybrid fitness and assortative mating 
strength and store the outcome of each simulation in a JLD2 file.

# Arguments
- `set_name::String`: the name assigned to the set of simulations.
- `intrinsic_R::Real`: the intrinsic growth rate.
- `K_total::Integer=40000`: the carrying capacity of the environment.
- `max_generations::Integer=1000`: the number of generations that the simulation will run.
- `total_loci::Integer=6`: the total number of loci in the genome.
- `female_mating_trait_loci=1:3`: the loci specifying the female's mate preference.
- `male_mating_trait_loci=1:3`: the loci specifying the male's mating trait.
- `hybrid_survival_loci=1:3`: the loci specifying the probability of survival to adulthood.
- `survival_fitness_method:String="epistasis"`: the method used to calculate the probability 
of survival to adulthood.
- `per_reject_cost=0`: the fitness loss of female per male rejected (due to search time, 
etc.). Can take values of 0 to 1.
- `sigma_disp=0.05`: the standard deviation of the normal distribution determining how far 
offspring will disperse from their mothers.
"""
function run_HZAM_set(
	set_name::String,
	total_loci::Int,
	female_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	male_mating_trait_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	hybrid_survival_loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	per_reject_cost::Real,
	w_hyb_set::Union{UnitRange{<:Real}, Vector{<:Real}},
	S_AM_set::Union{UnitRange{<:Real}, Vector{<:Real}};
	intrinsic_R::Real = 1.1, K_total::Int = 30000, max_generations::Int = 2000,
	survival_fitness_method::String = "epistasis", sigma_disp = 0.03,
)
	dir = mkpath(string(results_folder, "/", set_name))
	set_working_dir(dir)

	# Loop through the different simulation sets
	for i in eachindex(w_hyb_set)
		for j in eachindex(S_AM_set)
			w_hyb = w_hyb_set[i]
			S_AM = S_AM_set[j]

			run_name_half_gen =
				string(
					"HZAM_simulation_run_gen", trunc(max_generations / 2), "_SC",
					per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM,
				)

			println("Beginning: $run_name_half_gen")
			# run one simulation by calling the function defined above:
			outcome = run_one_HZAM_sim(
				w_hyb,
				S_AM,
				intrinsic_R;
				K_total,
				max_generations,
				total_loci,
				female_mating_trait_loci,
				male_mating_trait_loci,
				hybrid_survival_loci,
				survival_fitness_method,
				per_reject_cost,
				sigma_disp,
				do_plot = false,
				run_name = run_name_half_gen,
			)
			println("Ending: $run_name_half_gen")

			if !isnothing(outcome)
				run_name_full_gen =
					string(
						"HZAM_simulation_run_gen", max_generations, "_SC",
						per_reject_cost, "_Whyb", w_hyb, "_SAM", S_AM,
					)
				filename = string(dir, "/", run_name_full_gen, ".jld2")
				@save filename outcome
			end
		end # of S_AM loop
	end # of w_hyb loop   
end


"""
	set_results_folder(dir::String)

Set the output parent directory for all methods in summarize_data.jl.
"""
function set_results_folder(dir::String)
	global results_folder = mkpath(dir)
end

"""
	set_working_dir(dir::String)

Set the output directory for all methods in summarize_data.jl.
"""
function set_working_dir(dir::String)
	global working_dir = mkpath(dir)
end

"""
	plot_output_field(
		outcomes::Array{<:Real},
		sim_params::Array{DataAnalysis.SimParams}
	)

Create a heatmap of an output variable vs hybrid fitness and assortative mating.

# Arguments
- `outcomes::Array{:Real}`: the output from the simulation to be displayed.
- `sim_params::Array{<:DataAnalysis.SimParams}`: the simulation parameters resulting in the
 outcomes.
"""
function plot_output_field(
	outcomes::Array{<:Real},
	sim_params::Array{DataAnalysis.SimParams},
)
	output = [outcomes...]
	w_hybs = [[s.w_hyb for s in sim_params]...]
	S_AMs = [[s.S_AM for s in sim_params]...]

	w_hyb_set = sort(union(w_hybs))
	S_AM_set = sort(union(S_AMs))
	xticks = collect(1:length(S_AM_set))
	yticks = collect(1:length(w_hyb_set))

	xs = map(s -> indexin(s, S_AM_set)[1], S_AMs)
	ys = map(w -> indexin(w, w_hyb_set)[1], w_hybs)


	fontsize_theme = Theme(fontsize = 25)
	set_theme!(fontsize_theme)  # this sets the standard font size

	fig = Figure(resolution = (1800, 1200), figure_padding = 60, colormap = :grayC)

	ax = Axis(
		fig[1, 1],
		xlabel = "w_hyb",
		ylabel = "S_AM",
		xticks = (xticks, string.(S_AM_set)),
		yticks = (yticks, string.(w_hyb_set)),
		xticklabelsize = 20,
		yticklabelsize = 20,
		xlabelsize = 15,
		ylabelsize = 15,
	)
	hm = heatmap!(ax, xs, ys, output, colorrange = (0, 0.5))

	Colorbar(fig[:, 2], hm, ticklabelsize = 15)
	display(fig)
	readline()
	return fig
end

"""
	load_from_folder(dir::String)

Load the data from each simulation output file stored in the given directory into an array 
organized by hybrid fitness and assortative mating strength.
"""
function load_from_folder(dir::String)
	files = readdir(dir)
	outcome_array = Array{Union{DataAnalysis.OutputData, Missing}, 2}(
		undef, length(w_hyb_set), length(S_AM_set),
	)
	for i in eachindex(outcome_array)
		outcome_array[i] = missing
	end
	for file in files
		path = string(dir, "/", file)
		if occursin(".jld2", path) && occursin("gen1000", path)
			@load path outcome

			println(outcome.sim_params)

			if outcome.sim_params.S_AM in S_AM_set
				outcome_array[
					indexin(outcome.sim_params.w_hyb, w_hyb_set)[1],
					indexin(outcome.sim_params.S_AM, S_AM_set)[1],
				] = outcome
			end
		end
	end
	return outcome_array
end

"""
	load_from_folder(dir::String)

Load the data from each simulation output file stored in the given directory where hybrid 
fitness is one and store the output in a vector of tuples `(outcome, S_AM)`.
"""
function load_from_folder_whyb_1(dir::String)
	files = readdir(dir)
	outcome_array = Tuple{DataAnalysis.OutputData, <:Real}[]
	for file in files
		path = string(dir, "/", file)
		if occursin(".jld2", path) && occursin("gen1000", path) && occursin("Whyb1", path)
			if occursin(".jld2", path) && occursin("gen1000", path)
				alt_filepath = replace(path, "gen1000.0" => "gen2000")

				if isfile(alt_filepath)
					@load alt_filepath outcome
				else
					@load path outcome
				end

				println(outcome.sim_params)
				width = mean(outcome.hybrid_zone_width[end-19:end])
				if true #width < 20 * outcome.sim_params.sigma_disp
					push!(outcome_array, (outcome, outcome.sim_params.S_AM))
				end
			end
		end
	end
	return outcome_array
end


"""
	plot_fitnesses(fitnesses::Vector{<:Dict})

Produce a plot of fitnesses per phenotype over time.

# Arguments
- `fitnesses::Vector{<:Dict}`: number of offspring per phenotype each generation.
"""
function plot_fitnesses(fitnesses::Vector{<:Dict})
	dir = mkpath(string(results_folder, "/plots/"))
	@save string(dir, "fitness.jld2") fitnesses
	generations = collect(eachindex(fitnesses))
	num_generations = length(generations)
	phenotypes = collect(keys(fitnesses[1]))
	num_phenotypes = length(fitnesses[1])
	xs = vcat([fill(gen, num_phenotypes) for gen in eachindex(fitnesses)]...)
	ys = vcat([collect(keys(f)) for f in fitnesses]...)
	zs = vcat([collect(values(f)) for f in fitnesses]...)

	"""
	Determines which phenotype has the highest fitness each generation.
	"""
	function maximum_fitness()
		return map(
			f -> reduce((x, y) -> f[x] ≥ f[y] ? x : y, keys(f)),
			fitnesses,
		)
	end

	"""
	Compute the running average of the fitnesses for a given phenotype.
	"""
	function average_fitness(phenotype)
		n = 10
		average_fitnesses = []
		for i in eachindex(fitnesses)
			fitness_values =
				[fitnesses[j][phenotype] for j in i:min(length(fitnesses), i + n - 1)]
			push!(average_fitness, mean(fitness_values))
		end
		return average_fitnesses
	end

	fig = Figure(resolution = (1800, 1200), figure_padding = 60)
	ax = Axis(fig[1, 1], xlabel = "generation", ylabel = "phenotype")
	scatter!(ax, xs, ys)
	display(fig)
	display(heatmap(xs, ys, zs, colormap = :grayC))
	readline()
end

"""
	plot_population_tracking_data(filepath::String)

Create plots of the population size, hybrid zone width, hybrid index, and population overlap 
vs time.

# Arguments
- `filepath::String`:: the filepath for the file containing the outcome of the simulation.
"""
function plot_population_tracking_data(filepath::String)
	@load filepath outcome
	population_tracking_data = outcome.population_tracking_data
	sim_params = outcome.sim_params

	overlaps = [t.overlap for t in population_tracking_data]
	hybridity = [t.hybridity for t in population_tracking_data]
	widths = [t.width for t in population_tracking_data]
	populations = [t.population_size for t in population_tracking_data]

	xs = collect(eachindex(population_tracking_data))
	ax = Vector(undef, 4)
	points = Vector(undef, 4)

	fontsize_theme = Theme(fontsize = 40)
	set_theme!(fontsize_theme)  # set the standard font size
	fig = Figure(resolution = (1800, 1200), figure_padding = 60)
	# create the axis and labels
	ax[1] = Axis(
		fig[1, 1],
		ylabel = "# individuals",
		yticklabelsize = 30.0f0,
	)
	ax[2] = Axis(
		fig[2, 1],
		ylabel = "cline width",
		limits = ((nothing, nothing), (-0.05, 1.05)),
	)
	ax[3] = Axis(
		fig[3, 1],
		ylabel = "hybridity",
		yticklabelsize = 30.0f0,
	)
	ax[4] = Axis(
		fig[4, 1],
		ylabel = "overlap",
		xlabel = "generation",
		limits = ((nothing, nothing), (-0.05, 1.05)),
	)
	points[1] = scatter!(ax[1], xs, populations)
	points[2] = scatter!(ax[2], xs, widths)
	points[3] = scatter!(ax[3], xs, hybridity)
	points[4] = scatter!(ax[4], xs, overlaps)

	yspace = 1.5 * maximum(tight_yticklabel_spacing!, [ax[1], ax[2]])

	ax[1].yticklabelspace = yspace
	ax[2].yticklabelspace = yspace
	ax[3].yticklabelspace = yspace
	ax[4].yticklabelspace = yspace

	display(fig)

	readline()

	dir = mkpath(string(results_folder, "/plots/"))

	save(string(
		dir,
		"_S_AM",
		sim_params.S_AM,
		"_w_hyb",
		sim_params.w_hyb,
		".png",
		fig,
	))

	readline()
end

"""
	summarize_gene_correlations(source_dir::String; filename="gene_correlations")

Create plots of the gene correlations between different traits vs S_AM and w_hyb for 
different combinations of the same loci controlling different traits.

source_dir should point to a folder containing magic_preference.jld2, magic_cue.jld2, 
search_cost.jld2, and no_magic.jld2. Each file should contain an outcome array from a 
simulation set.

# Arguments
- `source_dir::String`: the folder where the outcomes of all the simulations are stored.
- `file_name="gene_correlations`: the name of the image file for the plot. 
"""
function summarize_gene_correlations(source_dir::String; filename = "gene_correlations")
	names = ["magic_preference", "magic_cue", "search_cost", "no_magic"]
	xlabelnames = ["Female mating trait", "Male mating trait", "Neutral trait"]
	ylabelnames = ["Magic preference", "Magic cue", "Search cost", "No pleiotropy"]
	xlabels = []
	ylabels = []
	correlations = Matrix(undef, 4, 3)

	magic_loci = [[1, 5], [3, 5], [5, 6], [5, 6]]
	fmt_loci = [[2], 1:2, 1:2, 1:2]
	mmt_loci = [3:4, [4], 3:4, 3:4]
	neutral_loci = [6:9, 6:9, 7:10, 7:10]

	for i in 1:4
		filename = string(source_dir, "/", names[i], ".jld2")
		@load filename sim_params outcome_array
		correlations[i, 1] =
			vcat(map(
				genotypes -> DataAnalysis.calc_trait_correlation(
					genotypes,
					magic_loci[i],
					fmt_loci[i],
				),
				outcome_array,
			)...)

		correlations[i, 2] =
			vcat(map(
				genotypes -> DataAnalysis.calc_trait_correlation(
					genotypes,
					magic_loci[i],
					mmt_loci[i],
				),
				outcome_array,
			)...)

		correlations[i, 3] =
			vcat(map(
				genotypes -> DataAnalysis.calc_trait_correlation(
					genotypes,
					magic_loci[i],
					neutral_loci[i],
				),
				outcome_array,
			)...)

	end

	w_hyb_set = string.([0.95, 0.9, 0.7, 0.5])
	S_AM_set = string.([1, 10, 100, 1000])
	xticks = [4, 3, 2, 1]
	yticks = [1, 2, 3, 4]

	xs = vcat(fill(xticks, 4)...)
	ys = vcat([fill(s, 4) for s in yticks]...)

	ax = Matrix(undef, 4, 4)
	hm = Matrix(undef, 4, 4)


	fontsize_theme = Theme(fontsize = 25)
	set_theme!(fontsize_theme)  # this sets the standard font size

	fig = Figure(
		resolution = (1800, 1200),
		figure_padding = 60,
		colorrange = (0, 1),
		colormap = :curl,
	)

	println(correlations)

	for j in 1:4
		for i in 1:3
			ax[j, i] = Axis(
				fig[j, i],
				xlabel = "w_hyb",
				ylabel = "S_AM",
				xticks = (xticks, w_hyb_set),
				yticks = (yticks, S_AM_set),
				xticklabelsize = 20,
				yticklabelsize = 20,
				xlabelsize = 15,
				ylabelsize = 15,
			)
			hm[j, i] = heatmap!(ax[j, i], xs, ys, correlations[j, i], colorrange = (-1, 1))
			println(correlations[j, i][16])
		end
		println("")
	end

	for i in 1:3
		push!(xlabels, Label(fig[0, i], xlabelnames[i], tellwidth = false))
	end

	for i in 1:4
		push!(
			ylabels,
			Label(fig[i, 0], ylabelnames[i], rotation = pi / 2, tellheight = false),
		)
	end

	Colorbar(fig[:, 4], hm[1, 1], label = "Pearson coefficient", ticklabelsize = 15)
	display(fig)
	dir = mkpath(string(results_folder, "/plots/"))
	save(string(dir, "/", filename, ".png"), fig)
	readline()
end
