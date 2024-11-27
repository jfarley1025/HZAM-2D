include("$(dirname(@__DIR__))/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using GLMakie # needed for making plots
using Statistics # needed for calculating standard deviation
GLMakie.activate!()

"The base of the figure."
f = Figure(size = (1000, 600))

"The figure layout."
layout = f[1, 1] = GridLayout()

"Titles for the subplots."
subplot_titles = ["A) 1 locus", "B) 3 loci"]

"Folders from which the simulation outcomes are loaded."
folders = [
	[
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/search_cost_outcomes/search_cost_cline_1locus",
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/search_cost_outcomes/search_cost_cline_sc_0.01_1locus",
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/search_cost_outcomes/search_cost_cline_sc_0.05_1locus",
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/search_cost_outcomes/search_cost_cline_sc_0.1_1locus",
	],
	[
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/Run2_2024July22/full_pleiotropy",
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/Run2_2024July22/low_reject_full_pleiotropy",
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/Run2_2024July22/high_reject_full_pleiotropy",
		"$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/search_cost_outcomes/search_cost_cline_sc_0.1_3loci",
	],
]

"Names of data series for the legend."
names = [
	"Search cost: 0",
	"Search cost: 0.01",
	"Search cost: 0.05",
	"Search cost: 0.1",
]

"Axes for the subplots."
ax = Vector{Makie.Axis}(undef, 2)

for j in 1:2
	# Set up the axis for the subplot
	ax[j] = Axis(
		layout[1, j],
		xscale = log10,
		xlabel = "Strength of assortative mating",
		ylabel = "Cline width",
		xlabelsize = 38,
		ylabelsize = 38,
		xticklabelsize = 35,
		yticklabelsize = 35,
		titlesize = 44,
		title = subplot_titles[j],
		aspect = 1,
		xgridvisible = false,
		ygridvisible = false,
		yticks = ([5, 10, 15, 20, 25], ["5σ", "10σ", "15σ", "20σ", "25σ"]),
		xticks = ([2, 10, 100, 1000], ["2×", "10×", "100×", "1000×"]),
	)

	# Add the cline width data for each search cost
	for i in 1:4
		outcome_array = HZAM.load_from_folder_whyb_1(folders[j][i])
		widths = [A.hybrid_zone_width[end-19:end] for (A, S_AM) in outcome_array]

		width = [mean(w) for w in widths] ./ Ref(0.03)
		S_AMs = [S_AM for (A, S_AM) in outcome_array]
		σs = [std(w) for w in widths] ./ Ref(0.03)
		perm = sortperm(S_AMs)

		lines!(ax[j], S_AMs[perm], width[perm], label = names[i], color = Makie.wong_colors()[i])
		errorbars!(ax[j], S_AMs[perm], width[perm], σs[perm], color = Makie.wong_colors()[i])
	end
end

ylims!(ax[2], 0.0, 30)

"Figure legend."
leg = Legend(layout[2, 1:2], ax[1], orientation = :horizontal, tellheight = true, labelsize = 38)

trim!(layout)
ax[1].tellwidth = true
linkyaxes!(ax[1], ax[2])
linkxaxes!(ax[1], ax[2])

display(f)
readline()
save("$(dirname(@__DIR__))/figures/search_cost_cline.png", f)