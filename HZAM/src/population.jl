"Types and functions used for managing population data (locations, genotypes, and growth rates)."
module Population
import ..DataAnalysis.calc_traits_additive
import QuadGK.quadgk

"The width of the (square) grid of zones that divides the population into more manageable 
chunks. Default is 10x10."
NUM_ZONES = 10

"""
    new_location(
        starting_range_x::Vector{Float32},
        starting_range_y::Vector{Float32}
    )::Tuple{Float32,Float32}

Randomly generate a location within the given range.

# Arguments 

- `starting_range_x::Vector{Float32}`: the x limits of the range
- `starting_range_y::Vector{Float32}`: the y limits of the range
"""
function new_location(
    starting_range_x::Vector{Float32},
    starting_range_y::Vector{Float32}
)::Tuple{Float32,Float32}
    x = rand() * (starting_range_x[2] - starting_range_x[1]) + starting_range_x[1]
    y = rand() * (starting_range_y[2] - starting_range_y[1]) + starting_range_y[1]
    return x, y
end

"""
    new_location(x::Float32, y::Float32, sigma_disp::Real)::Tuple{Float32,Float32}

Compute the location after dispersal of an offspring.

# Arguments 

- `x::Float32`: the x coordinate of the mother
- `y::Float32`: the y coordinate of the mother
- `sigma_disp::Real`: the standard deviation in the dispersal distance
"""
function new_location(x::Float32, y::Float32, sigma_disp::Real)::Tuple{Float32,Float32}
    new_x = -1
    new_y = -1
    # Keep generating new locations until there's one that's within the range.
    while ~(0 <= new_x < 1 && 0 <= new_y < 1) # Check if the location is within the range.
        dist = sigma_disp * randn()
        dir = rand() * 2 * pi
        new_x = x + dist * cos(dir)
        new_y = y + dist * sin(dir)
    end
    return new_x, new_y
end

"""
    Zone

A Zone stores the genotypes and locations for the individuals contained within 
a subset of the range.

# Fields
- `genotypes_F::Vector{Matrix{Int8}}`: the female genotypes. rows are alleles (row 1 from mother, row 2 from father) and columns are loci.
- `genotypes_M::Vector{Matrix{Int8}}`: the male genotypes.
- `x_locations_F::Vector{Float32}`: the x coordinates of the females
- `y_locations_F::Vector{Float32}`: the y coordinates of the females
- `x_locations_M::Vector{Float32}`: the x coordinates of the males
- `y_locations_M::Vector{Float32}`: the y coordinates of the males

# Constructors
```julia
- Zone(
    starting_N::Integer,
    total_loci::Integer,
    x_location::Float32,
    y_location::Float32,
    size::Float32,
    species::Real
)
- Zone(
    genotypes_F::Vector{Matrix{Int8}},
    genotypes_M::Vector{Matrix{Int8}},
    x_locations_F::Vector{Float32},
    y_locations_F::Vector{Float32},
    x_locations_M::Vector{Float32},
    y_locations_M::Vector{Float32}
)
```

# Details on behaviour of different constructors

The first constructor sets up the initial locations and genotypes of every 
individual in the zone.

The second constructor creates a new Zone using the genotypes and locations of the 
offspring.
"""
struct Zone
    genotypes_F::Vector{Matrix{Int8}}
    genotypes_M::Vector{Matrix{Int8}}
    x_locations_F::Vector{Float32}
    y_locations_F::Vector{Float32}
    x_locations_M::Vector{Float32}
    y_locations_M::Vector{Float32}

    function Zone(
        starting_N::Integer,
        total_loci::Integer,
        x_location::Float32,
        y_location::Float32,
        size::Float32,
        species::Real
    )
        species = Int8(species)

        N_half = trunc(Int, starting_N / 2) # the number of individuals in each sex

        zone_range_x = [x_location, min(x_location + size, 1.0f0)]
        zone_range_y = [y_location, min(y_location + size, 1.0f0)]
        # the geographic limits of the zone

        genotype = fill(species, 2, total_loci)

        genotypes = fill(genotype, N_half)

        x_locations_F = Vector{Float32}(undef, N_half)
        y_locations_F = Vector{Float32}(undef, N_half)
        x_locations_M = Vector{Float32}(undef, N_half)
        y_locations_M = Vector{Float32}(undef, N_half)

        for i in 1:N_half
            x_locations_F[i], y_locations_F[i] = new_location(zone_range_x, zone_range_y)
            x_locations_M[i], y_locations_M[i] = new_location(zone_range_x, zone_range_y)
        end

        new(
            genotypes,
            genotypes,
            x_locations_F,
            y_locations_F,
            x_locations_M,
            y_locations_M
        )
    end

    function Zone(
        genotypes_F::Vector{Matrix{Int8}},
        genotypes_M::Vector{Matrix{Int8}},
        x_locations_F::Vector{Float32},
        y_locations_F::Vector{Float32},
        x_locations_M::Vector{Float32},
        y_locations_M::Vector{Float32}
    )
        new(
            genotypes_F,
            genotypes_M,
            x_locations_F,
            y_locations_F,
            x_locations_M,
            y_locations_M
        )
    end
end

"""
    PopulationData

A PopulationData stores the genotypes, locations, and growth rates for all individuals in 
the simulation. The data is subdivided by 'zone' into a matrix for 
calculation efficiency.

# Fields
- `population::Matrix{Zone}`: all of the zones (containing the genotypes and locations) in the simulation.
- `growth_rates_F::Matrix{Vector{Float64}}`: the growth rates of every female for each zone.

# Constructors
```julia
- function PopulationData(
    K_total::Integer,
    total_loci::Integer,
    intrinsic_R::Real,
    sigma_comp::Real
)
- PopulationData(
    genotypes_daughters::Matrix{Vector{Matrix{Int8}}},
    genotypes_sons::Matrix{Vector{Matrix{Int8}}},
    x_locations_daughters::Matrix{Vector{Float32}},
    y_locations_daughters::Matrix{Vector{Float32}},
    x_locations_sons::Matrix{Vector{Float32}},
    y_locations_sons::Matrix{Vector{Float32}},
    K_total::Integer,
    sigma_comp::Real,
    intrinsic_R::Real
)
```
# Details on behaviour of different constructors

The first constructor sets up the initial population for the simulation.

The second constructor creates a new PopulationData using the genotypes and locations of the 
offspring.
"""
struct PopulationData
    population::Matrix{Zone}
    growth_rates_F::Matrix{Vector{Float64}}

    #=
    Compute the initial number of individuals in each zone.
    =#
    function calc_zone_population(K_total::Integer)
        floor(Int, K_total / (NUM_ZONES^2))
    end

    function PopulationData(
        K_total::Integer,
        total_loci::Integer,
        intrinsic_R::Real,
        sigma_comp::Real
    )
        # the number of individuals per zone (innitially constant throughout range)
        num_individuals_per_zone = calc_zone_population(K_total)

        # the locations of the zones along an axis
        spaced_locations = collect(0.0f0:Float32(1 / NUM_ZONES):0.99f0)

        zones = Matrix{Zone}(undef, NUM_ZONES, NUM_ZONES)

        # initialize the genotypes and locations for each zone
        for i in 1:NUM_ZONES
            for j in 1:NUM_ZONES
                zones[i, j] = Zone(
                    num_individuals_per_zone,
                    total_loci,
                    spaced_locations[i],
                    spaced_locations[j],
                    Float32(1 / NUM_ZONES),
                    spaced_locations[i] < 0.5 ? 0 : 1
                )
            end
        end

        # empty matrix for storing a list of growth rates of every female per zone
        growth_rates_F = Matrix{Vector{Float64}}(undef, NUM_ZONES, NUM_ZONES)

        # initialize the growth rates for every female
        for zone_index in CartesianIndices((1:NUM_ZONES, 1:NUM_ZONES))
            growth_rates_F[zone_index] = calc_growth_rates(zones,
                zone_index,
                K_total,
                sigma_comp,
                intrinsic_R)
        end

        new(zones, growth_rates_F)
    end

    function PopulationData(
        genotypes_daughters::Matrix{Vector{Matrix{Int8}}},
        genotypes_sons::Matrix{Vector{Matrix{Int8}}},
        x_locations_daughters::Matrix{Vector{Float32}},
        y_locations_daughters::Matrix{Vector{Float32}},
        x_locations_sons::Matrix{Vector{Float32}},
        y_locations_sons::Matrix{Vector{Float32}},
        K_total::Integer,
        sigma_comp::Real,
        intrinsic_R::Real
    )

        #= 
        Set up empty matrices for storing the zone data (genotypes and locations for every 
        individual in the zone) and the female growth rates.
        =#
        zones = Matrix{Zone}(undef, NUM_ZONES, NUM_ZONES)
        growth_rates_F = Matrix{Vector{Float64}}(undef, NUM_ZONES, NUM_ZONES)

        # initialize each zone with the offspring data
        for zone_index in CartesianIndices((1:NUM_ZONES, 1:NUM_ZONES))
            zones[zone_index] = Zone(
                genotypes_daughters[zone_index],
                genotypes_sons[zone_index],
                x_locations_daughters[zone_index],
                y_locations_daughters[zone_index],
                x_locations_sons[zone_index],
                y_locations_sons[zone_index]
            )
        end

        # calculate growth rates for every female in each zone
        for zone_index in CartesianIndices((1:NUM_ZONES, 1:NUM_ZONES))
            growth_rates_F[zone_index] = calc_growth_rates(zones,
                zone_index,
                K_total,
                sigma_comp,
                intrinsic_R)
        end

        new(zones, growth_rates_F)

    end
end

"""
    assign_zone(x::Float32, y::Float32)

Determine which zone a location falls in.
"""
function assign_zone(x::Float32, y::Float32)
    zone_x = convert(Integer, trunc(NUM_ZONES * x) + 1) # x index of the zone
    zone_y = convert(Integer, trunc(NUM_ZONES * y) + 1) # y index of the zone
    return CartesianIndex(
        min(zone_x, NUM_ZONES),
        min(zone_y, NUM_ZONES)
    )
end

"""
    calc_species_overlap(
        population::Matrix{Zone},
        max_dist::Real,
        sigma_comp::Real,
        loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
    )

Compute the total overlap between the two species.

# Arguments
- `population::Matrix{Zone}`: the genotypes and locations of all individuals.
- `max_dist::Real`: the distance cutoff for the density calculation.
- `sigma_comp::Real`: the standard deviation used for the density calculation.
- `loci::Union{UnitRange{<:Integer},Vector{<:Integer}}`: the loci used to determine species.
"""
function calc_species_overlap(
    population::Matrix{Zone},
    max_dist::Real,
    sigma_comp::Real,
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
    min_proportion = 0.1
    spaced_locations = collect(0.0f0:0.01f0:0.09f0)
    num_overlap = 0
    num_hybrid = 0
    num_A = 0
    num_B = 0

    function overlap_at_point(density_A, density_B, total_density)
        return density_A > total_density * min_proportion &&
               density_B > total_density * min_proportion
    end

    function hybrid_at_point(density_A, density_B, total_density)
        hybrid_density = total_density - density_A - density_B
        return hybrid_density > total_density * min_proportion
    end

    for i in 0:9
        for j in 0:9
            locations_x = Ref(Float32(i / 10)) .+ spaced_locations
            locations_y = Ref(Float32(j / 10)) .+ spaced_locations

            locations_x = vcat(fill(locations_x, 10)...)
            locations_y = vcat(fill.(locations_y, Ref(10))...)

            zone_index = CartesianIndex(i, j)
            density_A = calc_real_densities(
                population,
                zone_index,
                locations_x,
                locations_y,
                max_dist,
                sigma_comp;
                species=0,
                loci
            )

            density_B = calc_real_densities(
                population,
                zone_index,
                locations_x,
                locations_y,
                max_dist,
                sigma_comp;
                species=1,
                loci
            )

            total_density = calc_real_densities(
                population,
                zone_index,
                locations_x,
                locations_y,
                max_dist,
                sigma_comp
            )

            num_overlap += count(
                i -> overlap_at_point(density_A[i], density_B[i], total_density[i]),
                eachindex(locations_x)
            )

            num_hybrid += count(
                i -> hybrid_at_point(density_A[i], density_B[i], total_density[i]),
                eachindex(locations_x)
            )

            num_A += count(
                i -> density_A[i] > total_density[i]*min_proportion,
                eachindex(locations_x) 
            )

            num_B += count(
                i -> density_B[i] > total_density[i]*min_proportion,
                eachindex(locations_x) 
            )
        end
    end

    return num_overlap / 100^2, num_hybrid / 100^2, num_A / 100^2, num_B / 100^2
end


"""
    max_radius_squared(x::Float32, y::Float32, t::Real, max_dist::Real)

Compute the distance from a point along an angle to the limit of the range 
(defined by x∈[0,1), y∈[0,1)). 

The distance gets cut off at max_dist. Used in calculating the ideal densities.

# Arguments
- `x::Float32`: the starting x coordinate
- `y::Float32`: the starting y coordinate
- `t::Real`: the angle.
- `max_dist`: the maximum distance.
"""
function max_radius_squared(x::Float32, y::Float32, t::Real, max_dist::Real)
    if y > (1 - max_dist) && t < pi
        y_dist = (1 - y) / sin(t)
    elseif y < max_dist && t > pi
        y_dist = -y / sin(t)
    else
        y_dist = max_dist
    end

    if x > (1 - max_dist) && (t < pi / 2 || t > 3 * pi / 2)
        x_dist = (1 - x) / cos(t)
    elseif x < max_dist && (pi / 2 < t < 3 * pi / 2)
        x_dist = x / cos(pi - t)
    else
        x_dist = max_dist
    end

    return (min(x_dist, y_dist, max_dist))^2
end

"""
    calc_ideal_densities(
        K_total::Integer,
        sigma_comp::Real,
        x_locations_F::Vector{Float32},
        y_locations_F::Vector{Float32},
        max_dist::Real
    )

Compute expected local densities, based on geographically even distribution of individuals 
at carrying capacity.


Densities are calculated using the following integral:

``K_{total}\\int\\limits_0^{2\\pi}\\int\\limits_0^{0.03} 
\\exp(-\\frac{r^2}{2σ_{comp}^2})rdrdθ``

where the density function is equal to 0 outside of the range limits.

# Argument
- `K_total::Integer`: the carrying capacity
- `sigma_comp::Real`: the standard deviation for the normal curve used in calculating local density.
- `x_locations_F::Vector{Float32}`: the x coordinates where the ideal density is to be computed.
- `y_locations_F::Vector{Float32}`: the y coordinates where the ideal density is to be computed.
- `max_dist::Real`: the cutoff for the furthest away an individual can be and affect the density calculation. Should be 3x sigma_comp.
"""
function calc_ideal_densities(
    K_total::Integer,
    sigma_comp::Real,
    x_locations_F::Vector{Float32},
    y_locations_F::Vector{Float32},
    max_dist::Real
)
    function calc_ideal_density(x::Float32, y::Float32)
        if max_dist < x < (1 - max_dist) &&
           max_dist < y <= (1 - max_dist)
            # the ideal density is constant further than 3 sigma_comps from the range limits 
            return begin
                1 +
                K_total * 2 * pi *
                (1 - exp(-((max_dist^2) / 2) / (sigma_comp^2))) * (sigma_comp^2)
            end
        else
            # density function
            f(t) = exp(-(max_radius_squared(x, y, t, max_dist) / (2 * sigma_comp^2)))
            return begin
                1 +
                K_total * (sigma_comp^2) *
                (2 * pi -
                 quadgk(f, 0, 2 * pi)[1])
            end
        end
    end

    return calc_ideal_density.(x_locations_F, y_locations_F)
end

"""
    function calc_squared_distances(
        x_locations::Vector{Float32},
        y_locations::Vector{Float32},
        x::Float32,
        y::Float32,
        cutoff::Real
    )

Compute the squared distances from a set of locations to a single location. Discards 
distances above the given cutoff.

# Example

```jldoctest
julia> calc_squared_distances([0.4f0, 0f0, 0.5f0], [0.4f0, 0f0, 0.5f0],0.1)
2-element Vector{Float32}:
 0.0
 0.019999998
``` 
"""
function calc_squared_distances(
    x_locations::Vector{Float32},
    y_locations::Vector{Float32},
    x::Float32,
    y::Float32,
    cutoff::Real
)
    dif_x = x_locations .- Ref(x)
    dif_y = y_locations .- Ref(y)
    if cutoff == -1
        return dif_x .^ 2 .+ dif_y .^ 2
    else
        return filter!(x -> x < cutoff, dif_x .^ 2 .+ dif_y .^ 2)
    end
end

"""
    get_surrounding_zones(zone_index::CartesianIndex)

Return a list of the indices of the neighbouring zones and the given zone index.
"""
function get_surrounding_zones(zone_index::CartesianIndex)
    lower_left = max(zone_index - CartesianIndex(1, 1), CartesianIndex(1, 1))
    upper_right = min(
        zone_index + CartesianIndex(1, 1),
        CartesianIndex(NUM_ZONES, NUM_ZONES)
    )
    return lower_left:upper_right
end



"""
    calc_real_densities(
        population::Matrix{Zone},
        zone_index::CartesianIndex,
        x_locations_F::Vector{Float32},
        y_locations_F::Vector{Float32},
        max_dist::Real,
        sigma_comp::Real;
        species=-1,
        loci=[-1]
    )

Compute the population density at each female's location.

# Arguments
- `population::Matrix{Zone}`: the matrix of zones storing all the locations, and genotypes in the simulation.
- `zone_index::CartesianIndex`: the index of the focal zone.
- `x_locations_F::Vector{Float32}`: the x coordinates of every female in the focal zone.
- `y_locations_F::Vector{Float32}`: the y coordinates of every female in the focal zone.
- `max_dist::Real`: the distance cutoff for the density calculation.
- `sigma_comp::Real`: the standard deviation used in the density calculation.
"""
function calc_real_densities(
    population::Matrix{Zone},
    zone_index::CartesianIndex,
    x_locations_F::Vector{Float32},
    y_locations_F::Vector{Float32},
    max_dist::Real,
    sigma_comp::Real;
    species=-1,
    loci=[-1]
)
    neighbourhood = population[get_surrounding_zones(zone_index)]

    x_locations_all = vcat([[z.x_locations_F; z.x_locations_M] for z in neighbourhood]...)
    y_locations_all = vcat([[z.y_locations_F; z.y_locations_M] for z in neighbourhood]...)

    """
    Calculate the density at a given location
    """
    function calc_density(focal_x::Float32, focal_y::Float32)
        squared_distances = calc_squared_distances(
            x_locations_all,
            y_locations_all,
            focal_x,
            focal_y,
            max_dist^2
        )
        return sum(exp.(-squared_distances ./ Ref(2 * (sigma_comp^2))))
    end

    if species ≠ -1
        genotypes_all = vcat([[z.genotypes_F; z.genotypes_M] for z in neighbourhood]...)
        hybrid_indices = calc_traits_additive(genotypes_all, loci)
        filtered_indices = filter(i -> hybrid_indices[i] == species, eachindex(hybrid_indices))

        x_locations_all = x_locations_all[filtered_indices]
        y_locations_all = y_locations_all[filtered_indices]

        return calc_density.(x_locations_F, y_locations_F)
    end


    return calc_density.(x_locations_F, y_locations_F)
end

"""
    calc_growth_rates(
        population::Matrix{Zone},
        zone_index::CartesianIndex,
        K_total::Integer,
        sigma_comp::Real,
        intrinsic_R::Real
    )

Compute the female growth rates assuming no ecological difference.

# Arguments
- `population::Matrix{Zone}`: the matrix of zones storing all the locations, and genotypes in the simulation.
- `zone_index::CartesianIndex`: the index of the focal zone.
- `K_total::Integer`: the carrying capacity.
- `sigma_comp::Real`: the standard deviation for the normal curve used in calculating local density.
- `intrinsic_R::Real`: the intrinsic growth rate.
"""
function calc_growth_rates(
    population::Matrix{Zone},
    zone_index::CartesianIndex,
    K_total::Integer,
    sigma_comp::Real,
    intrinsic_R::Real
)
    # maximum distance that the density calculation takes into account
    max_dist = 3 * sigma_comp

    #= 
    Set up the expected local densities, based on a geographically even distribution of 
    individuals at the carrying capacity.
    =#
    ideal_densities = calc_ideal_densities(
        K_total,
        sigma_comp,
        population[zone_index].x_locations_F,
        population[zone_index].y_locations_F,
        max_dist
    )

    real_densities = calc_real_densities(
        population, zone_index,
        population[zone_index].x_locations_F,
        population[zone_index].y_locations_F,
        max_dist,
        sigma_comp
    )

    """
    Calculate the local growth rates according to discrete time logistic growth equation.
    """
    function logistic_growth_equation(ideal_densities, real_densities)
        return intrinsic_R .*
               ideal_densities ./ (ideal_densities .+ (real_densities .* (intrinsic_R - 1)))
    end

    return logistic_growth_equation(ideal_densities, real_densities)
end

"""
    genotype_mean(genotype::Matrix{<:Integer}, loci::Union{UnitRange{<:Integer},Vector{<:Integer}})

Compute the mean value of the genotype across the list of loci given.

# Example 
```jldoctest
julia> genotype_mean([0 1 0; 0 1 0], [1,3])
0.0
```
"""
function genotype_mean(
    genotype::Matrix{<:Integer},
    loci::Union{UnitRange{<:Integer},Vector{<:Integer}}
)
    return sum(genotype[:, loci]) / (2 * length(loci))
end

"""
    calc_survival_fitness_hetdisadvantage(genotype::Matrix{<:Integer}, w_hyb::Real)

Compute the survival fitness of an individual according to heterozygosity.

# Arguments
- `genotype:Matrix{<:Integer}`: the individual's genotype.
- `w_hyb::Real`: the simulation's hybrid fitness value.
"""
function calc_survival_fitness_hetdisadvantage(
    genotype::Matrix{<:Integer},
    w_hyb::Real
)::Float32
    num_loci = size(genotype, 2)
    s_per_locus = 1 - w_hyb^(1 / num_loci)  # loss in fitness due to each heterozygous locus 
    num_hetloci = sum(genotype[1, :] .≠ genotype[2, :])

    hetdisadvantage_fitness = (1 - s_per_locus)^num_hetloci
    return hetdisadvantage_fitness
end

"""
    calc_survival_fitness_epistasis(
        genotype::Matrix{<:Integer}, 
        hybrid_survival_loci::Union{UnitRange{<:Integer},Vector{<:Integer}}, 
        w_hyb::Real, 
        beta::Real=1
    )

Compute the survival fitness of an individual according to epistasis, with the beta 
parameter set to one as a default.

# Arguments
- `genotype:Matrix{<:Integer}`: the individual's genotype.
- `hybrid_survival_loci::Union{UnitRange{<:Integer},Vector{<:Integer}}`: the loci responsible for heterozygote disadvantage in computing the survival fitness.
- `w_hyb::Real`: the simulation's hybrid fitness value.
- `beta::Real`
"""
function calc_survival_fitness_epistasis(genotype::Matrix{<:Integer},
    hybrid_survival_loci::Union{UnitRange{<:Integer},Vector{<:Integer}},
    w_hyb::Real,
    beta::Real=1
)::Float32
    survival_HI = genotype_mean(genotype, hybrid_survival_loci)
    epistasis_fitnesses = 1 - (1 - w_hyb) * (4 * survival_HI * (1 - survival_HI))^beta
    return epistasis_fitnesses
end
end