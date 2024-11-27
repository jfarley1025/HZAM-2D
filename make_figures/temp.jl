include("$(dirname(@__DIR__))/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using CairoMakie
using Colors, ColorSchemes

"Colors associated with each phenotype (under 3 loci) for the viridis colormap."
colors = cgrad(:viridis, 7, categorical = true)

"Mean function."
mean(A) = sum(A) / length(A)

"""
	make_phenotype_subplot(g::Makie.GridLayout, filename::String, plotname::String)

Make the subplot showing how the phenotype frequencies change over time, the phenotype 
densities along a transect, and the location/phenotype of every individual at the end of the 
simulation.

# Arguments
- `g::Makie.GridLayout`: the part of the figure where the subplot is to be placed.
- `filename::String`: path to the simulation outcome.
- `plotname::String`: the title of the subplot.
"""
function make_phenotype_subplot(g::Makie.GridLayout, filename::String, plotname::String)
	colgap!(g, 10)
	rowgap!(g, 10)

	"Change the phenotypes array into a matrix."
	function reshape_output(phenotypes)
		output = Matrix(undef, length(phenotypes), length(phenotypes[1]))

		for i in eachindex(phenotypes)
			for j in eachindex(phenotypes[1])
				output[i, j] = phenotypes[i][j]
			end
		end
		output
	end

	@load filename outcome

    fmt_phenotype_counts, mmt_phenotype_counts = outcome.population_tracking_data

	"The phenotype frequency over time and transect sub-subplots are grid layouts since each 
	row of the plot is a separate heatmap because the colormaps need to be different."
	g1 = g[1, 1] = GridLayout()
	g2 = g[2, 1] = GridLayout()
	g3 = g[1, 2] = GridLayout()
	g4 = g[2, 2] = GridLayout()

	"Axis showing the preference phenotype frequencies over time."
	ax1 = Axis(
		g[1, 1],
		ylabel = "Preference",
		yticks = ([1, 7], ["0", "1"]),
		ylabelsize = 22,
		xticklabelsize = 20,
		yticklabelsize = 20,
		yminorticksvisible = true,
		yminorticks = [4],
	)

	"Axis showing the cue phenotype frequencies over time."
	ax2 = Axis(
		g[2, 1],
		xlabel = "generation",
		ylabel = "Cue",
		yticks = ([1, 7], ["0", "1"]),
		xlabelsize = 22,
		ylabelsize = 22,
		xticklabelsize = 20,
		yticklabelsize = 20,
		yminorticksvisible = true,
		yminorticks = [4],
	)

	"Axis showing the preference densities over a transect at the end of the simulation."
	ax3 = Axis(
		g[1, 2],
		yminorticksvisible = true,
		xminorticksvisible = true,
		yminorticks = [4],
		xminorticks = [51],
	)


	"Axis showing the cue densities over a transect at the end of the simulation."
	ax4 = Axis(
		g[2, 2],
		xlabel = "x",
		xlabelsize = 22,
		xticklabelsize = 20,
		yticklabelsize = 20,
		yminorticksvisible = true,
		xminorticksvisible = true,
		yminorticks = [4],
		xminorticks = [51],
	)

	"Axis showing the locations and preferences of each individual in the simulation."
	ax5 = Axis(
		g[1, 3],
		aspect = 1,
	)

	"Axis showing the locations and cues of each individual in the simulation."
	ax6 = Axis(
		g[2, 3],
		aspect = 1,
	)

	pd = outcome.population_data
	genotypes = [
		vcat([d.genotypes_F for d in pd.population]...)
		vcat([d.genotypes_M for d in pd.population]...)
	]
	x_locations = [
		vcat([d.x_locations_F for d in pd.population]...)
		vcat([d.x_locations_M for d in pd.population]...)
	]
	y_locations = [
		vcat([d.y_locations_F for d in pd.population]...)
		vcat([d.y_locations_M for d in pd.population]...)
	]


	"Add the preference phenotype densities to ax3."
	plot_phenotype_distribution(
		ax3,
		g3,
		genotypes,
		outcome.sim_params.female_mating_trait_loci,
		x_locations,
		y_locations,
	)

	"Add the cue phenotype densities to ax4."
	plot_phenotype_distribution(
		ax4,
		g4,
		genotypes,
		outcome.sim_params.male_mating_trait_loci,
		x_locations,
		y_locations,
	)

	"Add the locations and preferences of each individual to ax5."
	plot_population(
		ax5,
		genotypes,
		x_locations,
		y_locations,
		outcome.sim_params.female_mating_trait_loci,
	)

	"Add the locations and cues of each individual to ax5."
	plot_population(
		ax6,
		genotypes,
		x_locations,
		y_locations,
		outcome.sim_params.male_mating_trait_loci,
	)


	hidexdecorations!(ax1, ticks = false, minorticks = false)
	hideydecorations!(ax4, ticks = false, minorticks = false)
	hidedecorations!(ax3, ticks = false, minorticks = false)
	hidedecorations!(ax5, ticks = true)
	hidedecorations!(ax6, ticks = true)

	ax2.xticks = [1000, 2000]
	linkxaxes!(ax1, ax2)
	linkyaxes!(ax1, ax2)
	ax1.xticks = [1000, 2000]

	"Add the preference phenotype frequencies over time to ax1."
    println(length(fmt_phenotype_counts))
	heatmap!(ax1, reshape_output(fmt_phenotype_counts), colormap = Reverse(:grayC))
	composite_heatmap(g1, reshape_output(fmt_phenotype_counts))

	"Add the cue phenotype frequencies over time to ax1."
	heatmap!(ax2, reshape_output(mmt_phenotype_counts), colormap = Reverse(:grayC))
	composite_heatmap(g2, reshape_output(mmt_phenotype_counts))

	"Add the subplot title."
	Label(g[1, :, Top()], plotname, valign = :bottom,
		font = :bold,
		padding = (0, 0, 5, 0),
		fontsize = 22)

	colsize!(g, 1, Auto(2))
	colsize!(g, 2, Auto(2))
	colsize!(g, 3, Auto(1))
	trim!(g)

	g
end

"""
	plot_population(
		ax ::Makie.Axis,
		genotypes::Vector{<:Matrix{<:Real}},
		x_locations::Vector{<:Real},
		y_locations::Vector{<:Real},
		loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	)

Add the locations and phenotype of each individual in the simulation to a sub-subplot.

# Arguments
-`ax ::Makie.Axis`: axis where the points are to be added.
-`genotypes::Vector{<:Matrix{<:Real}}`: genotypes of every individual in the simulation.
-`x_locations::Vector{<:Real}`: x coordinate of every individual in the simulation.
-`y_locations::Vector{<:Real}`: y coordinate of every individual in the simulation.
-`loci::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: loci used to calculate trait values.
"""
function plot_population(
	ax::Makie.Axis,
	genotypes::Vector{<:Matrix{<:Real}},
	x_locations::Vector{<:Real},
	y_locations::Vector{<:Real},
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
)
	xlims!(ax, (0, 1))
	ylims!(ax, (0, 1))
	hybrid_indices = HZAM.DataAnalysis.calc_traits_additive(genotypes, loci)

	scatter!(ax, x_locations, y_locations, color = hybrid_indices, markersize = 1,
		colorrange = (0, 1), colormap = :viridis)

	ax.xticks = ([0, 1], ["0", "1"])
	ax.yticks = ([0, 1], ["0", "1"])
end

"""
	plot_phenotype_distribution(
		ax::Makie.Axis,
		g::Makie.GridLayout,
		genotypes::Vector{<:Matrix{<:Integer}},
		loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		x_locations::Vector{<:Real},
		y_locations::Vector{<:Real},
	)

Add the phenotype densities along the y=0.5 transect to a GridLayout. The GridLayout where 
the coloured heatmap is added is placed over a black and white heatmap to keep the 
formatting. The GridLayout is needed because the coloured heatmap is really a composite of 
7 heatmaps (one for each row/phenotype/colour).

# Arguments
-`ax ::Makie.Axis`: axis over which the heatmap is to be added.
-`g::Makie.GridLayout`: the GridLayout onto which the heatmap is to be added.
-`genotypes::Vector{<:Matrix{<:Real}}`: genotypes of every individual in the simulation.
-`loci::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: loci used to calculate trait values.
-`x_locations::Vector{<:Real}`: x coordinate of every individual in the simulation.
-`y_locations::Vector{<:Real}`: y coordinate of every individual in the simulation.
"""
function plot_phenotype_distribution(
	ax::Makie.Axis,
	g::Makie.GridLayout,
	genotypes::Vector{<:Matrix{<:Integer}},
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	x_locations::Vector{<:Real},
	y_locations::Vector{<:Real},
)
	output = Array{Float64, 2}(undef, 2 * length(loci) + 1, 101)
	total_densities = Array{Float64}(undef, 101)
	for x in 0:1:100
		dif_x = x_locations .- Ref(x / 100)
		dif_y = y_locations .- Ref(0.5)
		squared_distances = dif_x .^ 2 .+ dif_y .^ 2
		d = sum(exp.(-squared_distances ./ Ref(2 * (0.01^2))))
		total_densities[x+1] = d
	end

	genotypes = [g[:, loci] for g in genotypes]
	for phenotype in 0:(2*length(loci))
		indices = findall(g -> count(x -> x == 1, g) == phenotype, genotypes)
		xs = x_locations[indices]
		ys = y_locations[indices]
		gs = genotypes[indices]
		for x in 0:1:100
			dif_x = xs .- Ref(x / 100)
			dif_y = ys .- Ref(0.5)
			squared_distances = dif_x .^ 2 .+ dif_y .^ 2
			d = sum(exp.(-squared_distances ./ Ref(2 * (0.01^2)))) / total_densities[x+1]
			output[phenotype+1, x+1] = d
		end
	end
	heatmap!(
		ax,
		output',
		colormap = Reverse(:grayC),
		colorrange = (0, 1),
	)
	composite_heatmap(g, output')
	ax.xticks = ([1, 101], ["0", "1"])

	ax.yticks = ([1, 2 * length(loci) + 1], ["0", "1"])
end

"""
	compare_phenotype_plots()

Create a plot showing four example simulation outcomes.
"""
function compare_phenotype_plots()
	f = Figure(size = (600, 300))

	ga = f[1, 1] = GridLayout()

	path = "matching.jld2"
	make_phenotype_subplot(ga, path, "Matching")

	display(f)
	save("$(dirname(@__DIR__))/figures/temp_phenotypes.png", f)
end

"""
	composite_heatmap(g::Makie.GridLayout, output)

Create a coloured heatmap from a matrix where each row uses a different colormap.

# Arguments
- `g::Makie.GridLayout`: the scene on which the coloured heatmap is to lie.
- `output`: the data used for the heatmap.
"""
function composite_heatmap(g::Makie.GridLayout, output)
	ax = Vector(undef, 7)
	for i in 1:7
		ax[i] = Axis(g[8-i, 1])
		m = reshape(output[:, i], length(output[:, i]), 1)
		heatmap!(ax[i], m, colormap = cgrad([:white, colors[i]]), colorrange = (0, 1))
		hidedecorations!(ax[i])
		hidespines!(ax[i])
	end
	rowgap!(g, 0)
	colgap!(g, 0)
	trim!(g)
end

compare_phenotype_plots()
