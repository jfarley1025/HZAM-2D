@testset "calc_sigmoid_curve" begin
	locations_x = collect(0:0.01:0.99)
	hybrid_indices = [fill(0, 50); fill(1, 50)]

	s, width = DataAnalysis.calc_sigmoid_curve(locations_x, hybrid_indices)
	midpoint = DataAnalysis.spaced_locations[argmin(abs.(s .- Ref(0.5)))]
	@test length(s) == 101
	@test midpoint ≈ 0.5
	@test s[1] == 0
	@test s[end] == 1
	@test width < 0.05
end

@testset "sigmoid" begin
	xs = collect(0:0.01:1)

	s = DataAnalysis.sigmoid(xs, [0, 1])

	@test s[1] ≈ 0.5
	@test s[101] ≈ 1 / (1 + exp(-1))
end


@testset "calc_traits_additive" begin
	genotypes = [[0 0 0; 0 0 0], [1 1 1; 1 1 1]]

	@test DataAnalysis.calc_traits_additive(genotypes, 1:3) == [0, 1]

	genotypes = [[1 0 0; 0 0 0], [0 1 0; 1 0 0], [1 1 1; 1 1 1]]

	@test DataAnalysis.calc_traits_additive(genotypes, 1:2) == [0.25, 0.5, 1]
end


@testset "calc_cline_width" begin
	pd1 = Population.PopulationData(10000, 2, 1.1, 0.01)
	cline_width_1 = DataAnalysis.calc_cline_width(pd1, 1:2, 0.1:0.2:0.9)
	@test cline_width_1 < 0.05

	pd2 = Population.PopulationData(0, 2, 1.1, 0.01)
	cline_width_2 = DataAnalysis.calc_cline_width(pd2, 1:2, 0.1:0.2:0.9)

	@test cline_width_2 > 1
end

@testset "calc_transect" begin
	genotypes1 = vcat(fill([0 0; 0 0], 5000), fill([1 1; 1 1], 5000))
	genotypes2 = fill([0 0; 1 1], 10000)
	x_locations = collect(0.0001:0.0001:1)
	y_locations = rand(10000)

	avg_hi_1, sigmoid_1, cline_width_1 = DataAnalysis.calc_transect(
		genotypes1,
		x_locations,
		y_locations,
		1:2,
		0.01,
		0.5,
	)
	@test avg_hi_1[30] < 0.01
	@test 1 - avg_hi_1[60] < 0.01
	@test sigmoid_1[30] < 0.01
	@test 1 - sigmoid_1[60] < 0.01
	@test cline_width_1 < 0.05

	avg_hi_2, sigmoid_2, cline_width_2 = DataAnalysis.calc_transect(
		genotypes2,
		x_locations,
		y_locations,
		1:2,
		0.01,
		0.5,
	)
	@test all(avg_hi_2 .== 0.5)
	@test all((sigmoid_2 .- 0.5) .< 0.01)
	@test cline_width_2 > 1
end

@testset "calc_transects" begin
	pd1 = Population.PopulationData(10000, 2, 1.1, 0.01)
	hi_avgs, sigmoid_curves, widths = DataAnalysis.calc_transects(pd1, 1:2, 0.1:0.2:0.9)
	@test all(widths .< 0.05)
    for s in sigmoid_curves
        @test s[30] < 0.01
        @test 1 - s[60] < 0.01
    end
    for h in hi_avgs
        @test h[30] < 0.01
        @test 1 - h[60] < 0.01
    end

	pd2 = Population.PopulationData(0, 2, 1.1, 0.01)
	hi_avgs, sigmoid_curves, widths = DataAnalysis.calc_transects(pd2, 1:2, 0.1:0.2:0.9)

	for h in hi_avgs
        @test all(h .== 0)
    end

    for w in widths
        @test w > 1
    end
end