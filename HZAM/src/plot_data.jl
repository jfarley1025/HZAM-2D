"Functions for plotting data while the simulation is running."
module PlotData

export create_population_plot, update_population_plot, create_gene_plot, update_gene_plot


#gr()  # use GR backend for graphs
using Colors, ColorSchemes
import ColorSchemes.plasma
#using Plots.PlotMeasures  # needed for plot margin adjustment
#using GLMakie
using ..DataAnalysis

# Activate the GLMakie backend
#GLMakie.activate!()

#GLMakie.activate!(inline=false) # set up the plot to display in its own window

# The figure that updates every generation.
global fig

# Axes on which the locations are plotted.
global ax

# Points showing the location and hybrid index of every individual.
global points

"""
    function create_population_plot(
        hybrid_indices_functional::Vector{<:Real},
        x_locations::Vector{Float32},
        y_locations::Vector{Float32},
        save_plot::Bool
    )

Initialize the plot of locations and hybrid indices.

# Arguments
- `hybrid_indices_functional::Vector{<:Real}`: list of the hybrid index (value between 0 and 1) of 
every individual.
- `x_locations::Vector{Float32}`: the x coordinate of every individual. 
- `y_locations::Vector{Float32}`: the y coordinate of every individual.
- `save_plot::Bool`: true if the plot is to be saved to a PNG file.
"""
function create_population_plot(
    hybrid_indices_functional::Vector{<:Real},
    x_locations::Vector{Float32},
    y_locations::Vector{Float32},
    save_plot::Bool
)
    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # set the standard font size
    global fig = Figure(resolution=(1800, 1200), figure_padding=60)
    # create the axis and labels
    global ax = Axis(
        fig[1, 1],
        xlabel="location_x",
        ylabel="location_y",
        title=string("HZAM simulation, generation = ", 0),
        xticklabelsize=45,
        yticklabelsize=45,
        titlegap=30
    )
    # set the limits for the plotted area
    xlims!(-0.03, 1.03)
    ylims!(-0.03, 1.03)

    # add the location of every individual to the plot
    global points = scatter!(
        ax,
        x_locations,
        y_locations,
        color=hybrid_indices_functional,
        markersize=10
    )
    if save_plot
        dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse_no_pleiotropy")
        save(string(dir, "/", 1, ".png"), fig)
    end

    display(fig)

    return fig
end

"""
    update_population_plot(
        hybrid_indices_functional::Vector,
        locations::Vector,
        generation::Integer,
        save_plot::Bool
    )

Update the existing plot with new locations and hybrid indices.

# Arguments
- `hybrid_indices_functional::Vector`: list of the hybrid index (value between 0 and 1) of 
every individual.
- `locations::Vector`: the location of every individual (must be in the same order as the 
hybrid indices).
- `generation::Integer`: the number of elapsed generations.
- `save_plot::Bool`: true if the plot is to be saved to a PNG file.
"""
function update_population_plot(
    hybrid_indices_functional::Vector,
    x_locations::Vector,
    y_locations::Vector,
    generation::Integer,
    save_plot::Bool
)
    delete!(ax, points) # remove the old points from the plot

    # add the location of every individual to the plot
    global points = scatter!(
        ax,
        x_locations,
        y_locations,
        color=hybrid_indices_functional,
        markersize=10
    )

    ax.title = string("HZAM simulation, generation = ", generation)


    println("generation: ", generation, "; individuals: ", length(x_locations))
    println("")
    
    if save_plot
        dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse_no_pleiotropy")
        save(string(dir, "/", generation, ".png"), fig)
    end
    display(fig)
end

"""
    function count_genotypes(
        range::Union{UnitRange{<:Integer},Vector{<:Integer}},
        genotypes::Vector{<:Matrix{<:Integer}}
    )

Count the number of individuals with each genotype for a given loci range.
"""
function count_genotypes(
    range::Union{UnitRange{<:Integer},Vector{<:Integer}},
    genotypes::Vector{<:Matrix{<:Integer}}
)
    num_loci = length(range)
    indices_range = (0:(1/(2*num_loci)):1)
    hybrid_indices = DataAnalysis.calc_traits_additive(genotypes, range)
    return map(i -> count(x -> x ≈ i, hybrid_indices), indices_range)
end

"""
    get_expected(x::Real, n::Integer)

Compute the expected number of individuals with a given trait hybrid index if the trait is 
controlled by 4 loci. 

# Arguments
- `x::Real`: the hybrid index for the trait of interest.
- `n::Integer`: the total population size.
"""
function get_expected(x::Real, n::Integer)
    return (n * binomial(8, Int(8 * x))) / 2^8
end

"""
    create_gene_plot(
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::NamedTuple,
        save_plot::Bool
    )

Initialize the plot of the phenotype frequencies for each loci range.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes of every individual.
- `loci::NamedTuple`: the loci range for each trait.
- `save_plot::Bool`: true if the plot is to be saved to a file.
"""
function create_gene_plot(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::NamedTuple,
    save_plot::Bool
)
    println(loci)

    plot_titles = ["Neutral trait", "Female mating trait", "Male mating trait",
        "Hybrid survival trait", "Expected"]
    global ax = Vector(undef, 5)
    global points = Vector(undef, 5)

    fontsize_theme = Theme(fontsize=60)
    set_theme!(fontsize_theme)  # set the standard font size
    global fig = Figure(resolution=(1800, 1200), figure_padding=60)

    # initialize the axis for each subplot
    for i in 1:5
        global ax[i] = Axis(
            fig[Int(ceil(i / 3)), (i-1)%3+1],
            xlabel="trait value",
            ylabel="# individuals",
            title=plot_titles[i],
            yticklabelsize=20,
            xticklabelsize=20
        )
    end

    hybrid_survival_indices = DataAnalysis.calc_traits_additive(genotypes, loci.hybrid_survival)

    species_A_indices = findall(h -> h == 0, hybrid_survival_indices)
    species_B_indices = findall(h -> h == 1, hybrid_survival_indices)
    other_indices = setdiff(eachindex(genotypes), union(species_A_indices, species_B_indices))

    species_A_output = count_genotypes.(values(loci), Ref(genotypes[species_A_indices]))
    species_B_output = count_genotypes.(values(loci), Ref(genotypes[species_B_indices]))
    other_output = count_genotypes.(values(loci), Ref(genotypes[other_indices]))

    for i in 3:6
        num_loci = length(values(loci)[i])
        range = (0:(1/(2*num_loci)):1)
        xs = [range; range; range]
        stk = [fill(1, length(species_A_output[i])); fill(3, length(species_B_output[i])); fill(2, length(other_output[i]))]
        ys = [species_A_output[i]; species_B_output[i]; other_output[i]]

        global points[i-2] = barplot!(
            ax[i-2],
            xs,
            ys,
            color=stk,
            stack=stk,
            colorrange=(1.1, 2.9),
            highclip=:orange,
            lowclip=:blue)
    end

    global points[5] = barplot!(
        ax[5],
        (0:0.125:1),
        get_expected.((0:0.125:1), Ref(length(genotypes))),
        color=:blue
    )
    println("generation 1")
    display(fig)

    if save_plot
        dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse_no_pleiotropy_high_search_cost")
        save(string(dir, "/", 1, ".png"), fig)
    end
end

"""
    update_gene_plot(
        genotypes::Vector{<:Matrix{<:Integer}},
        loci::NamedTuple,
        generation::Integer,
        save_plot::Bool
    )

Update the plot of the phenotype frequencies for each loci range.

# Arguments
- `genotypes::Vector{<:Matrix{<:Integer}}`: the genotypes of every individual.
- `loci::NamedTuple`: the loci range for each trait.
- `generation::Integer`: the generation # the simulation is on.
- `save_plot::Bool`: true if the plot is to be saved to a file.
"""
function update_gene_plot(
    genotypes::Vector{<:Matrix{<:Integer}},
    loci::NamedTuple,
    generation::Integer,
    save_plot::Bool
)

    [delete!(ax[i], points[i]) for i in 1:5] # remove the old points from the plot

    hybrid_survival_indices = DataAnalysis.calc_traits_additive(genotypes, loci.hybrid_survival)

    species_A_indices = findall(h -> h == 0, hybrid_survival_indices)
    species_B_indices = findall(h -> h == 1, hybrid_survival_indices)
    other_indices = setdiff(eachindex(genotypes), union(species_A_indices, species_B_indices))

    species_A_output = count_genotypes.(values(loci), Ref(genotypes[species_A_indices]))
    species_B_output = count_genotypes.(values(loci), Ref(genotypes[species_B_indices]))
    other_output = count_genotypes.(values(loci), Ref(genotypes[other_indices]))

    for i in 3:6
        num_loci = length(values(loci)[i])
        range = (0:(1/(2*num_loci)):1)
        xs = [range; range; range]
        stk = [fill(1, length(species_A_output[i])); fill(3, length(species_B_output[i])); fill(2, length(other_output[i]))]
        ys = [species_A_output[i]; species_B_output[i]; other_output[i]]

        global points[i-2] = barplot!(
            ax[i-2],
            xs,
            ys,
            color=stk,
            stack=stk,
            colorrange=(1.1, 2.9),
            highclip=:orange,
            lowclip=:blue)
    end

    global points[5] = barplot!(
        ax[5],
        (0:0.125:1),
        get_expected.((0:0.125:1), Ref(length(genotypes))),
        color=:blue
    )
    println(string("generation: ", generation))

    if save_plot
        dir = mkpath("HZAM_Sym_Julia_results_GitIgnore/plots/gene_timelapse_no_pleiotropy_high_search_cost")
        save(string(dir, "/", generation, ".png"), fig)
    end
end
end # end of PlotData module