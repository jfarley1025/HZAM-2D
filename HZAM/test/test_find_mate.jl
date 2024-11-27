using Random

@testset "choose_closest_male_no_zones" begin
    elig_M = collect(1:1000)
    female_x, female_y = 0.5f0, 0.5f0

    x_locations_M = Float32.(rand(1000))
    y_locations_M = Float32.(rand(1000))
    male_x, male_y = 0.50001f0, 0.49999f0
    x_locations_M[579] = male_x
    y_locations_M[579] = male_y

    closest_male =
        Mating.choose_closest_male_from_zone(
            elig_M,
            x_locations_M,
            y_locations_M,
            female_x,
            female_y
        )

    @test x_locations_M[closest_male] == male_x
    @test y_locations_M[closest_male] == male_y

    for i in 1:100
        male_x, male_y = 0.55f0, 0.45f0
        x_locations_M = Float32.(rand(1000) .* Ref(0.4))
        y_locations_M = Float32.(rand(1000) .* Ref(0.4) .+ Ref(0.55))

        i = min(Int(round(rand() * 1000)), 1)

        x_locations_M[i] = male_x
        y_locations_M[i] = male_y

        closest_male =
            Mating.choose_closest_male_from_zone(
                elig_M,
                x_locations_M,
                y_locations_M,
                female_x,
                female_y
            )

        @test x_locations_M[closest_male] == male_x
        @test y_locations_M[closest_male] == male_y
    end
end

@testset "choose_closest_male_one_zone" begin
    # locations of the zones along an axis
    intervals = collect(0.0f0:0.1f0:0.99f0)

    zones = Matrix{Population.Zone}(undef, 10, 10)
    for i in 1:10
        for j in 1:10
            zones[i, j] = Population.Zone(
                1000,
                10,
                intervals[i],
                intervals[j],
                Float32(1 / 10),
                0
            )
        end
    end

    x_location_mother = 0.55f0
    y_location_mother = 0.55f0

    elig_M = Dict{CartesianIndex,Vector{Int64}}()

    for i in eachindex(IndexCartesian(), zones)
        elig_M[i] = 1:length(zones[i].x_locations_M)
    end
    female_zone_index = Population.assign_zone(x_location_mother, y_location_mother)

    closest_male, zone_index =
        Mating.choose_closest_male(
            zones, 
            [female_zone_index],
            elig_M, x_location_mother, y_location_mother, 0.1f0
        )

    distance = Mating.distance(
        zones[zone_index].x_locations_M[closest_male],
        zones[zone_index].y_locations_M[closest_male],
        x_location_mother,
        y_location_mother
    )

    @test distance < 0.01

    @test zone_index == female_zone_index
end

@testset "choose_closest_male_multiple_zones" begin
    # locations of the zones along an axis
    intervals = collect(0.0f0:0.1f0:0.99f0)
    zones = Matrix{Population.Zone}(undef, 10, 10)

    for i in 1:10

        for i in 1:10
            for j in 1:10
                zones[i, j] = Population.Zone(
                    1000,
                    10,
                    intervals[i],
                    intervals[j],
                    Float32(1 / 10),
                    0
                )
            end
        end

        x_location_mother = 0.45f0
        y_location_mother = 0.45f0

        x_locations_M = vcat(zones[i, i].x_locations_M, x_location_mother)
        y_locations_M = vcat(zones[i, i].x_locations_M, y_location_mother)

        genotypes_F = zones[i, i].genotypes_F
        genotypes_M = zones[i, i].genotypes_M
        push!(genotypes_M, genotypes_M[1])
        x_locations_F = zones[i, i].x_locations_F
        y_locations_F = zones[i, i].y_locations_F

        zones[i, i] = Population.Zone(
            genotypes_F,
            genotypes_M,
            x_locations_F,
            y_locations_F,
            x_locations_M,
            y_locations_M
        )

        elig_M = Dict{CartesianIndex,Vector{Int64}}()

        for i in eachindex(IndexCartesian(), zones)
            elig_M[i] = 1:length(zones[i].x_locations_M)
        end
        neighbourhood =
            filter(e -> length(elig_M[e]) > 0, CartesianIndex(1, 1):CartesianIndex(10, 10))

        closest_male, zone_index =
            Mating.choose_closest_male(
                zones,
                neighbourhood,
                elig_M,
                x_location_mother,
                y_location_mother,
                1.0f0
            )

        @test zones[zone_index].x_locations_M[closest_male] == x_location_mother
        @test zones[zone_index].y_locations_M[closest_male] == y_location_mother
        @test zone_index == CartesianIndex(i, i)
    end
end

@testset "choose_closest_male_empty_zone" begin
    zone = Population.Zone(0, 0, 0.0f0, 0.0f0, 0.1f0, 0)
    zones = [zone zone; zone zone]
    elig_M = Dict{CartesianIndex,Vector{Int64}}()
    elig_M[CartesianIndex(1, 1)] = []

    closest_male = Mating.choose_closest_male(
        zones,
        [CartesianIndex(1, 1)],
        elig_M,
        0.5f0,
        0.6f0,
        0.5f0
    )
    @test closest_male == (-1, -1)
end



@testset "generate_offspring_genotype" begin
    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 0 1; 1 0 1]]

    @test Mating.generate_offspring_genotype(genotypes[1], genotypes[2]) == [0 0 0; 1 1 1]
    @test Mating.generate_offspring_genotype(genotypes[1], genotypes[1]) == [0 0 0; 0 0 0]
    @test Mating.generate_offspring_genotype(genotypes[3], genotypes[1]) == [1 0 1; 0 0 0]
end


@testset "calc_match_strength_S_AM1" begin
    female_mating_trait_loci = 2:3
    male_mating_trait_loci = 1:2
    pref_SD = sqrt(-1 / (2 * log(1 / (1 + 10^(-15)))))

    pref_SD = Inf

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 0 0 0]]

    match_strength11 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[1],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength13 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength23 =
        Mating.calc_match_strength(
            genotypes[2],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength12 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[2],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength33 =
        Mating.calc_match_strength(
            genotypes[3],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    @test match_strength11 == 1
    @test match_strength33 == 1
    @test match_strength23 ≈ 1
    @test match_strength12 ≈ 1
    @test match_strength13 ≈ 1
end

@testset "calc_match_strength_S_AM_inf" begin
    female_mating_trait_loci = 2:3
    male_mating_trait_loci = 1:2

    pref_SD = 10^(-15)

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 0 0 0]]

    match_strength11 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[1],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength13 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength23 =
        Mating.calc_match_strength(
            genotypes[2],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength12 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[2],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength33 =
        Mating.calc_match_strength(
            genotypes[3],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    @test match_strength11 == 1
    @test match_strength33 == 1
    @test match_strength23 ≈ 0
    @test match_strength12 ≈ 0
    @test match_strength13 ≈ 0
end

@testset "calc_match_strength_S_AM_10" begin
    female_mating_trait_loci = 2:3
    male_mating_trait_loci = 1:2

    pref_SD = sqrt(-1 / (2 * log(1 / 10)))

    genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1], [1 1 1; 0 0 0]]

    match_strength11 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[1],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength13 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength23 =
        Mating.calc_match_strength(
            genotypes[2],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength12 =
        Mating.calc_match_strength(
            genotypes[1],
            genotypes[2],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    match_strength33 =
        Mating.calc_match_strength(
            genotypes[3],
            genotypes[3],
            pref_SD,
            female_mating_trait_loci,
            male_mating_trait_loci
        )

    @test match_strength11 == 1
    @test match_strength33 == 1
    @test match_strength23 ≈ 0.5623413251903492
    @test match_strength12 ≈ 0.1
    @test match_strength13 ≈ 0.5623413251903492
end

@testset "distance" begin
    @test Mating.distance(0.5f0, 0.5f0, 0.5f0, 0.5f0) == 0
    @test Mating.distance(0.5f0, 0.5f0, 0.0f0, 0.0f0) ≈ sqrt(0.5)
    @test Mating.distance(0.5f0, 0.5f0, 1.0f0, 1.0f0) ≈ sqrt(0.5)
    @test Mating.distance(0.0f0, 0.0f0, 1.0f0, 1.0f0) ≈ sqrt(2)
end