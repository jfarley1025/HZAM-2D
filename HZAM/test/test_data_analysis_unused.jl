@testset "calc_all_trait_correlations" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 1], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)

    for tc in trait_correlations
        @test tc[3] == 1
    end

    genotypes = [[1 1 1; 1 1 1], [0 1 0; 0 1 0], [0 1 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)
    expected = [1, 0.9333, 0.93333, 0.9333, 1, 1, 0.9333, 1, 1]

    trait_correlations = [tc[3] for tc in trait_correlations]

    @test all(abs.(trait_correlations .- expected) .< 0.001)
end

@testset "calc_trait_correlation" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    @test abs(0.6708 - DataAnalysis.calc_trait_correlation(genotypes, [1, 3], [2])) < 0.0001
    @test abs(1 - DataAnalysis.calc_trait_correlation(genotypes, 1:3, 1:3)) < 0.0001
    @test abs(0.4472 - DataAnalysis.calc_trait_correlation(genotypes, [1], [3])) < 0.0001
    @test abs(0.8407 - DataAnalysis.calc_trait_correlation(genotypes, 1:2, 2:3)) < 0.0001

    genotypes = [[1 0 1; 0 1 1], [0 1 1; 1 0 1], [0 0 1; 1 1 1]]

    correlation = DataAnalysis.calc_trait_correlation(genotypes, [1, 2, 3], [3])

    @test correlation == 0
end

@testset "calc_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 1)

    @test linkage_diseq == 1

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 2)

    @test linkage_diseq ≈ 1 / 9

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 3)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 3, 1)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 2, 3)

    @test linkage_diseq ≈ 1 / 2
end

@testset "calc_average_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]
    loci = (overall=1:3, l1=1:2, l2=[1, 3])

    average_linkage_diseq = DataAnalysis.calc_average_linkage_diseq(genotypes, loci)

    expected = [2576 / 5184, 615 / 1296, 103 / 216, 205 / 432, 5 / 9, 125 / 288, 103 / 216, 125 / 288, 9 / 16]

    averages = vcat([a[3] for a in average_linkage_diseq]...)

    @test averages ≈ expected
end

@testset "calc_length" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [[fill(0, i); fill(1, 100 - i)] for i in 40:5:85]
    sigmoid_curves = DataAnalysis.calc_sigmoid_curve.(Ref(locations_x), hybrid_indices)

    length = DataAnalysis.calc_length(sigmoid_curves)

    @test abs(length - 1.106) < 0.001

    hybrid_indices = fill([fill(0, 50); fill(1, 50)], 10)
    sigmoid_curves = DataAnalysis.calc_sigmoid_curve.(Ref(locations_x), hybrid_indices)

    length = DataAnalysis.calc_length(sigmoid_curves)

    @test length ≈ 1
end

@testset "calc_position" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [[fill(0, i); fill(1, 100 - i)] for i in 40:5:85]
    sigmoid_curves = DataAnalysis.calc_sigmoid_curve.(Ref(locations_x), hybrid_indices)

    position = DataAnalysis.calc_position(sigmoid_curves)

    @test position ≈ 0.621
end

@testset "calc_variance" begin
    positions = collect(0:0.01:0.99)

    @test DataAnalysis.calc_variance(positions) ≈ 0.01

    positions = zeros(500)

    @test DataAnalysis.calc_variance(positions) == 0
end

@testset "count_phenotypes_at_loci" begin
    genotypes = [fill([0 1; 0 1], 100); fill([1 0; 0 0], 100)]

    c1 = DataAnalysis.count_phenotypes_at_loci(genotypes, [1])
    c2 = DataAnalysis.count_phenotypes_at_loci(genotypes, [1, 2])

    @test c1[0.0f0] == 100
    @test c1[0.5f0] == 100
    @test c1[1.0f0] == 0
    @test c2[0.0f0] == 0
    @test c2[0.25f0] == 100
    @test c2[0.5f0] == 100
    @test c2[0.75f0] == 0
    @test c2[1.0f0] == 0
end

@testset "find_fixed_alleles" begin
    genotypes = [[1 1 1; 1 1 1], [1 1 1; 1 1 1], [1 0 0; 1 0 0], [1 0 1; 1 1 0]]

    @test DataAnalysis.find_fixed_alleles(genotypes) == [1]

    genotypes = fill([1 1 1; 1 1 1], 50)

    @test DataAnalysis.find_fixed_alleles(genotypes) == [1, 2, 3]

    genotypes = [[1 1; 1 1], [0 0; 0 0]]

    @test DataAnalysis.find_fixed_alleles(genotypes) == []
end

@testset "calc_chi_squared" begin
    genotypes = [[1 1 1; 1 1 1], [1 0 1; 1 1 1], [1 1 0; 1 1 0], [0 0 0; 0 0 0],
        [0 1 0; 1 0 1]]

    x = DataAnalysis.calc_chi_squared(genotypes, 1)

    @test abs(x - 1.3719) < 0.001

    x = DataAnalysis.calc_chi_squared(genotypes, 2)

    @test abs(x - 0.1389) < 0.001

    x = DataAnalysis.calc_chi_squared(genotypes, 3)

    @test abs(x - 1.8) < 0.001
end


@testset "calc_all_trait_correlations" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 1], [0 0 0; 0 0 0], [0 0 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)

    for tc in trait_correlations
        @test tc[3] == 1
    end

    genotypes = [[1 1 1; 1 1 1], [0 1 0; 0 1 0], [0 1 0; 0 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [1 1 1; 1 1 1]]

    trait_correlations = DataAnalysis.calc_all_trait_correlations(genotypes, loci)
    expected = [1, 0.9333, 0.93333, 0.9333, 1, 1, 0.9333, 1, 1]

    trait_correlations = [tc[3] for tc in trait_correlations]

    @test all(abs.(trait_correlations .- expected) .< 0.001)
end

@testset "calc_trait_correlation" begin
    loci = NamedTuple{(:l1, :l2, :l3),NTuple{3,Any}}((1:2, 3:3, [1, 3]))
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    @test abs(0.6708 - DataAnalysis.calc_trait_correlation(genotypes, [1, 3], [2])) < 0.0001
    @test abs(1 - DataAnalysis.calc_trait_correlation(genotypes, 1:3, 1:3)) < 0.0001
    @test abs(0.4472 - DataAnalysis.calc_trait_correlation(genotypes, [1], [3])) < 0.0001
    @test abs(0.8407 - DataAnalysis.calc_trait_correlation(genotypes, 1:2, 2:3)) < 0.0001

    genotypes = [[1 0 1; 0 1 1], [0 1 1; 1 0 1], [0 0 1; 1 1 1]]

    correlation = DataAnalysis.calc_trait_correlation(genotypes, [1, 2, 3], [3])

    @test correlation == 0
end

@testset "calc_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 1)

    @test linkage_diseq == 1

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 2)

    @test linkage_diseq ≈ 1 / 9

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 1, 3)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 3, 1)

    @test linkage_diseq ≈ 1 / 8

    linkage_diseq = DataAnalysis.calc_linkage_diseq(genotypes, 2, 3)

    @test linkage_diseq ≈ 1 / 2
end

@testset "calc_average_linkage_diseq" begin
    genotypes = [[1 1 1; 1 1 0], [0 0 0; 0 0 0], [1 0 0; 1 0 0], [0 0 0; 0 0 0],
        [1 1 1; 1 1 1], [0 1 0; 0 1 1]]
    loci = (overall=1:3, l1=1:2, l2=[1, 3])

    average_linkage_diseq = DataAnalysis.calc_average_linkage_diseq(genotypes, loci)

    expected = [2576 / 5184, 615 / 1296, 103 / 216, 205 / 432, 5 / 9, 125 / 288, 103 / 216, 125 / 288, 9 / 16]

    averages = vcat([a[3] for a in average_linkage_diseq]...)

    @test averages ≈ expected
end









@testset "calc_width" begin
    locations_x = collect(0.001:0.01:0.999)
    hybrid_indices = [fill(0, 50); fill(1, 50)]
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    width = DataAnalysis.calc_width(sigmoid_curve)

    @test abs(width) <= 0.005

    hybrid_indices = vcat(fill([0, 1], 50)...)
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    width = DataAnalysis.calc_width(sigmoid_curve)

    @test width ≈ 1


    hybrid_indices = collect(0:0.01:0.99)
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    width = DataAnalysis.calc_width(sigmoid_curve)

    @test abs(width - 0.9) <= 0.01
end



@testset "calc_sigmoid_curves" begin
    locations_x = collect(0.0f0:0.01f0:0.99f0)
    hybrid_indices = [fill(0, 50); fill(1, 50)]

    locations_x = vcat(fill(locations_x, 20)...)
    locations_y = vcat([fill(Float32((i - 1) * 0.05) + 0.001f0, 100) for i in 1:20]...)

    hybrid_indices = vcat(fill(hybrid_indices, 20)...)

    @test length(locations_x) == length(hybrid_indices)

    sigmoid_curves = DataAnalysis.calc_sigmoid_curves(locations_x, locations_y, hybrid_indices)
    s = sigmoid_curves[1]

    midpoint = DataAnalysis.spaced_locations[argmin(abs.(s .- 0.5))]
    @test maximum(sigmoid_curves[1].-sigmoid_curves[10]) <0.01
    @test length(s) == 1001
    @test midpoint ≈ 0.495
    @test s[1] == 0
    @test s[1000] == 1
end

@testset "PopulationTrackingData" begin
    genotypes = [[0 0; 1 1], [0 1; 0 0], [0 0; 0 0], [1 1; 1 1]]
    locations_x = [0.1f0, 0.6f0, 0.2f0, 0.7f0]
    locations_y = fill(0.51f0, 4)

    loci = NamedTuple{(:overall, :functional, :neutral, :female_mating_trait, :male_mating_trait, :competition_trait, :hybrid_survival),NTuple{7,Any}}(([1, 2, 3,], [1, 2], [3], 1:2, 1:2, 1:2, 1:2))

    tracking_data = DataAnalysis.PopulationTrackingData(genotypes, locations_x, locations_y, (1:2), 0)

    @test tracking_data.population_size == 4

    @test tracking_data.hybridity ≈ 1.5 / 4

    @test tracking_data.overlap ≈ 0

    @test tracking_data.width > 0.95
end

@testset "sort_locations" begin
    A = [0, 0, 0.5, 1, 0.76, 0.23]
    bin_size = 0.5

    sorted_indices = DataAnalysis.sort_locations(A, bin_size)

    @test sorted_indices == [[1, 2, 6], [3, 5]]

    bin_size = 0.1

    sorted_indices = DataAnalysis.sort_locations(A, bin_size)

    @test sorted_indices == [[1, 2], [], [6], [], [], [3], [], [5], [], []]
end

@testset "sort_y" begin
    @test DataAnalysis.sort_y([0.01, 0.5, 0.24, 0.9]) ==
          [[1], [], [], [], [3], [], [], [], [], [], [2], [], [], [], [], [], [], [], [4], []]
end
@testset "calc_bimodality_in_range" begin
    locations_x = collect(0:0.01:0.99)
    hybrid_indices = [fill(0, 50); fill(1, 50)]

    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)

    bimodality = DataAnalysis.calc_bimodality_in_range(
        sigmoid_curve,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality == 1

    hybrid_indices = [fill(0, 49); [0.5]; fill(1, 50)]

    bimodality = DataAnalysis.calc_bimodality_in_range(
        sigmoid_curve,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality ≈ 4 / 5
end

@testset "calc_bimodality_overall" begin
    locations_x = collect(0:0.01:0.99)
    hybrid_indices = [fill(0, 50); fill(1, 50)]
    sorted_indices = [collect(((i-1)*100+1):(((i-1)*100+1)+99)) for i in 1:10]
    sigmoid_curve = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)
    sigmoid_curves = fill(sigmoid_curve, 10)

    locations_x = vcat(fill(locations_x, 10)...)
    hybrid_indices = vcat(fill(hybrid_indices, 10)...)

    bimodality = DataAnalysis.calc_bimodality_overall(
        sigmoid_curves,
        sorted_indices,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality == 1

    h0 = [fill(0, 50); fill(1, 50)]
    h1 = [fill(0, 49); [0.95]; fill(1, 50)]

    hybrid_indices = vcat(h0, h1, h0, h1, h0, h1, h0, h1, h0, h1)

    bimodality = DataAnalysis.calc_bimodality_overall(
        sigmoid_curves,
        sorted_indices,
        locations_x,
        hybrid_indices,
        0.05
    )

    @test bimodality ≈ 0.9
end


@testset "average_width" begin
    locations_x = collect(0.001:0.01:0.999)
    h0 = [fill(0, 50); fill(1, 50)]
    h1 = vcat(fill([0, 1], 50)...)
    s0 = DataAnalysis.calc_sigmoid_curve(locations_x, h0)
    s1 = DataAnalysis.calc_sigmoid_curve(locations_x, h1)

    avg = DataAnalysis.average_width(vcat(fill([s0, s1], 5)...))

    @test abs(0.5 - avg) < 0.005
end