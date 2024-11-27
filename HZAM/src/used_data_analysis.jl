"Parameters for fitting sigmoid curves to the clines."
global initial_par = [0.0, 1.0]

"Evenly spaced locations across the range (one-dimensional)."
global spaced_locations = collect(Float32, 0:0.01:1)

"""
	SimParams

A SimParams stores the parameters used in a simulation.

# Fields
- `intrinsic_R::Real`: the intrinsic growth rate.
- `w_hyb::Real`: the hybrid fitness.
- `S_AM::Real`: the strength of assortative mating.
- `K_total::Integer`: the carrying capacity of the environment.
- `max_generations::Integer`: the number of generations the simulation ran for.
- `sigma_disp::Real`: the standard deviation for dispersal distance.
- `total_loci::Integer`: the number of loci in the genotypes.
- `female_mating_trait_loci`: the loci controlling mating preference.
- `male_mating_trait_loci`: the loci controlling mating cue.
- `hybrid_survival_loci`: the loci controlling hybrid fitness.
- `per_reject_cost::Real`: the loss in fitness incurred by a female after rejecting one male.
"""
struct SimParams
	intrinsic_R::Real
	w_hyb::Real
	S_AM::Real
	K_total::Integer
	max_generations::Integer
	sigma_disp::Real
	total_loci::Integer
	female_mating_trait_loci::Any
	male_mating_trait_loci::Any
	hybrid_survival_loci::Any
	per_reject_cost::Real
end


"""
	OutputData

An OutputData stores all of the key data from a simulation run.

# Fields
- `sim_params::SimParams`: the parameters of the simulation.
- `population_data`: the genotypes and locations for all individuals in the simulation.
- `hybrid_zone_width::Vector{<:Real}`: the cline width associated with the male mating trait.
- `population_overlap::Vector{<:Real}`: the proportion of the range occupied by males of both mating trait phenotypes.
- `bimodality::Real`: the extent to which phenotypically pure individuals occur in the hybrid zone
- `population_tracking_data::Vector{PopulationTrackingData}`: the population size, hybridity, overlap, and cline width over time.
- `overlap::Real`: the proportion of the range where both species occur.
"""
struct OutputData
	sim_params::SimParams
	population_data::Any
	hybrid_zone_width::Union{Real, Vector{<:Real}}
	population_overlap::Union{Real, Vector{<:Real}}
	bimodality::Real
	population_tracking_data::Any
end

"""
	calc_sigmoid_curve(locations_x::Vector{<:Real}, hybrid_indices::Vector{<:Real})

Compute a sigmoid curve that models hybrid index vs location on the x axis.

Return the output of the sigmoid curve function on 1000 evenly spaced locations.
"""
function calc_sigmoid_curve(locations_x::Vector{<:Real}, hybrid_indices::Vector{<:Real})
	# calculate the parameters for the curve fitting
	fit = curve_fit(sigmoid, locations_x, hybrid_indices, initial_par)
	# compute the values of the sigmoid curve at evenly spaced x values
	return sigmoid(spaced_locations, fit.param), 1 / (fit.param[2] / 4)
end


"""
	sigmoid(x::Vector{<:Real}, p::Vector{<:Real})

Compute the y value for a sigmoid given the x value, centre, and maximum slope.

# Arguments
- `x::Real`: the value along the x axis where the sigmoid is to be calculated.
- `p::Vector{<:Real}`: the fit parameters where p[1] is the cline centre and p[2] is 
the maximum slope.
"""
function sigmoid(x::Vector{<:Real}, p::Vector{<:Real})
	1 ./ (1 .+ exp.(-p[2] .* (x .- p[1])))
end

"""
	calc_traits_additive(
		genotypes::Vector{<:Matrix{<:Integer}},
		loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
	)::Vector{Float32}

Compute the mean values of the genotypes passed to it at the given loci. Used to determine 
trait values in an additive way.
"""
function calc_traits_additive(
	genotypes::Vector{<:Matrix{<:Integer}},
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
)::Vector{Float32}

	N = length(genotypes)
	# calculate the mean value of the genotype across the list of loci given
	function mean(genotype, loci)
		return sum(genotype[:, loci]) / (2 * length(loci))
	end

	traits = map(x -> mean(genotypes[x], loci), 1:N)
	return traits
end

"""
	calc_transect(
		genotypes::Vector{<:Matrix{<:Integer}}, 
		x_locations::Vector{<:Real}, 
		y_locations::Vector{<:Real}, 
		loci::Union{UnitRange{<:Integer},Vector{<:Integer}}, 
		sigma_comp::Real, 
		y_coord::Real
	)

Compute the cline width along a horizontal transect.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}` : the genotypes of every individual in the sim.
- `x_locations::Vector{<:Real}` : the x coordinate of every individual in the sim.
- `y_locations::Vector{<:Real}` : the y coordinate of every individual in the sim.
- `loci::Union{UnitRange{<:Integer},Vector{<:Integer}}` : the focal loci.
- `sigma_comp::Real` : the standard deviation used for the density calculation.
- `y_coord::Real` : the y coordinate of the transect along which cline width is measured.
"""
function calc_transect(
	genotypes::Vector{<:Matrix{<:Integer}},
	x_locations::Vector{<:Real},
	y_locations::Vector{<:Real},
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	sigma_comp::Real,
	y_coord::Real,
)
	spaced_locations = collect(Float32, 0:0.01:1)
	hybrid_indices = calc_traits_additive(genotypes, loci)
	function average_hybrid_index(x, y)
		dif_x = x_locations .- Ref(x)
		dif_y = y_locations .- Ref(y)

		squared_distances = dif_x .^ 2 .+ dif_y .^ 2
		weights = exp.(-squared_distances ./ Ref(2 * (sigma_comp^2)))
		return sum(weights) == 0 ? 0 : sum(hybrid_indices .* weights) / sum(weights)
	end

	avg_hybrid_indices = average_hybrid_index.(spaced_locations, Ref(y_coord))
	sigmoid, width = calc_sigmoid_curve(spaced_locations, avg_hybrid_indices)
	return avg_hybrid_indices, sigmoid, width
end

"""
	calc_transects(
		pd,
		loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		y_coords::Union{StepRangeLen{<:Real}, Vector{<:Real}},
	)

Compute a sigmoidal cline model for the specified loci at horizontal transects with the 
specified y coordinates. Returns the weighted average hybrid index at evenly spaced 
locations along the transect, the sigmoidal curves modeling the transects, and the cline 
widths.

# Arguments
- `pd`: the PopulationData containing the genotypes, x coordinates, and y coordinates of 
all individuals.
- `y_coords`: the y coordinates of the transects.
- `loci::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: the loci specifying the trait of 
interest.
"""
function calc_transects(
	pd,
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	y_coords::Union{StepRangeLen{<:Real}, Vector{<:Real}},
)
	hybrid_index_avgs = []
	sigmoid_curves = []
	widths = []

	sigma_comp = 0.01

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

	for y in y_coords
		avg_hybrid_indices, sigmoid, width = calc_transect(genotypes, x_locations, y_locations, loci, sigma_comp, y)
		push!(hybrid_index_avgs, avg_hybrid_indices)
		push!(sigmoid_curves, sigmoid)
		push!(widths, width)
	end

	return hybrid_index_avgs, sigmoid_curves, widths
end


"""
	calc_cline_width(
		pd,
		loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
		y_coords::Union{StepRangeLen{<:Real}, Vector{<:Real}},
	)

Compute the cline width by averaging the width along transects at the given y coordinates.

# Arguments
- `pd`: the PopulationData containing the genotypes, x coordinates, and y coordinates of 
all individuals.
- `y_coords`: the y coordinates of the transects used to compute the width.
- `loci::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: the loci specifying the trait of 
interest.
"""
function calc_cline_width(
	pd,
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	y_coords::Union{StepRangeLen{<:Real}, Vector{<:Real}},
)
	hybrid_index_avgs, sigmoid_curves, widths = calc_transects(pd, loci, y_coords)
	return mean(widths)
end

"""
	calc_bimodality_on_transect(
		pd,
		sigmoid_curve::Vector{<:Real},
		y_coord::Real,
		loci::Union{UnitRange{<:Integer}, Vector{<:Integer}}
	)

Compute the bimodality along a transect by locating the midpoint and calculating the 
density of phenotypically pure individuals there.

# Arguments
- `pd`: the PopulationData containing the genotypes, x coordinates, and y coordinates of 
all individuals.
- `sigmoid_curve::Vector{<:Real}`: a sigmoid curve fit to the cline along the transect.
- `y_coord::Real`: the y coordinate of the transect.
- `loci::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: the loci specifying the trait of 
interest.
"""
function calc_bimodality_on_transect(
	pd,
	sigmoid_curve::Vector{<:Real},
	y_coord::Real,
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
)
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

	hybrid_indices = calc_traits_additive(genotypes, loci)

	y_coord = Float32(y_coord)
	midpoint = Float32(spaced_locations[argmin(abs.(sigmoid_curve .- Ref(0.5)))])

	"Compute the phenotype density at a point. `species` specifies which phenotype is of 
	interest and if `species=-1 then every phenotype is included."
	function calc_density(x, y; species = -1)
		if species != -1
			indices = findall(h -> h == species, hybrid_indices)
		else
			indices = eachindex(hybrid_indices)
		end

		dif_x = x_locations[indices] .- Ref(x)
		dif_y = y_locations[indices] .- Ref(y)
		squared_distances = dif_x .^ 2 .+ dif_y .^ 2

		return sum(exp.(-squared_distances ./ Ref(2 * (0.01^2))))
	end

	density = calc_density(
		midpoint,
		y_coord,
	)[1]

	density_A = calc_density(
		midpoint,
		y_coord,
		species = 0,
	)[1]

	density_B = calc_density(
		midpoint,
		y_coord,
		species = 1,
	)[1]

	return (density_A + density_B) / density

end

"""
	calc_bimodality(
		pd,
		loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
	)

Compute the bimodality of the hybrid zone by calculating the bimodality along 9 transects.

# Arguments
- `pd`: the PopulationData containing the genotypes, x coordinates, and y coordinates of 
all individuals.
- `loci::Union{UnitRange{<:Integer}, Vector{<:Integer}}`: the loci specifying the trait of 
interest.
"""
function calc_bimodality(
	pd,
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
)
	hybrid_index_avgs, sigmoid_curves, widths = calc_transects(pd, loci, 0.1:0.1:0.9)
	bimodalities =
		calc_bimodality_on_transect.(
			Ref(pd),
			sigmoid_curves,
			collect(0.1:0.1:0.9),
			Ref(loci),
		)

	return mean(bimodalities)
end

"""
	function count_phenotypes_at_loci(
		genotypes::Vector{<:Matrix{<:Integer}},
		loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
	)

Return the number of individuals with each phenotype for the given loci range.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes of every individual
- `loci::Union{UnitRange{<:Integer},Vector{<:Integer}}`: the loci range over which the phenotypes are calculated.
"""
function count_phenotypes_at_loci(
	genotypes::Vector{<:Matrix{<:Integer}},
	loci::Union{UnitRange{<:Integer}, Vector{<:Integer}},
)
	hybrid_indices = calc_traits_additive(genotypes, loci)
	n = 2 * length(loci)

	output = []

	for i in 0:n
		push!(output, count(x -> x ≈ i / n, hybrid_indices) / length(hybrid_indices))
	end

	return output
end