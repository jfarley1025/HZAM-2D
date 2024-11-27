include("$(dirname(@__DIR__))/HZAM/src/HZAM.jl")
import .HZAM
using JLD2 # needed for saving / loading data in Julia format
using CairoMakie # used for plotting
using Colors

@load "$(dirname(@__DIR__))/HZAM-2D_Julia_results_GitIgnore/example_outcomes/magic_preference" outcome

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

hybrid_indices = HZAM.DataAnalysis.calc_traits_additive(genotypes, outcome.sim_params.male_mating_trait_loci)

fig = Figure(resolution = (1600, 1000))
# create the axis and labels
ax = Axis(
	fig[1, 1],
	xlabel = "x",
	ylabel = "y",
	xlabelsize = 30,
	ylabelsize = 30,
	xticklabelsize = 30,
	yticklabelsize = 30,
)

g = fig[1, 2] = GridLayout()
# set the limits for the plotted area
xlims!(-0.03, 1.03)
ylims!(-0.03, 1.03)

# add the location of every individual to the plot
points = scatter!(
	ax,
	x_locations,
	y_locations,
	color = hybrid_indices,
	markersize = 8,
)


Colorbar(fig[:, 0], points, ticklabelsize = 30, flipaxis = false, label = "Male mating trait", labelsize = 30)

# y coordinates of each transect
y_coords = collect(0.1:0.2:0.9)

# calculate the hybrid indices and fit curves along each transect
hybrid_index_avgs, sigmoid_curves, widths = HZAM.DataAnalysis.calc_transects(
	pd,
	outcome.sim_params.male_mating_trait_loci,
	y_coords,
)

sigmoid_axis = []

# add plot for each of the 5 transects
for i in 1:5
	push!(
		sigmoid_axis,
		Axis(
			g[6-i, 1],
			xlabel = "x",
			xlabelsize = 30,
			xticklabelsize = 30,
			xgridvisible = false,
			ygridvisible = false,
		),
	)

	scatter!(
		sigmoid_axis[i],
		HZAM.spaced_locations,
		hybrid_index_avgs[i],
		markersize = 8,
	)
	lines!(
		sigmoid_axis[i],
		HZAM.spaced_locations,
		sigmoid_curves[i],
		linewidth = 10,
		color = (:gray, 0.5))
	text!(
		Point.(0.1, 0.8),
		text = "y = $(y_coords[i])",
		align = (:center, :center),
		color = :black,
		fontsize = 30,
	)
	hideydecorations!(sigmoid_axis[i])
	if i != 1
		hidexdecorations!(sigmoid_axis[i])
	end
end
rowgap!(g, 0)
colsize!(fig.layout, 2, Relative(3 / 7))

display(fig)

save(string("$(dirname(@__DIR__))/figures/methods_fig.png"), fig)
