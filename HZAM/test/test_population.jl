@testset "initialize_population_1" begin
	n = Population.NUM_ZONES
	K_total = 2 * (n^2)
	ecolDiff = 0
	total_loci = 3
	intrinsic_R = 1.1
	sigma_comp = 0.01

	pd = Population.PopulationData(K_total,
		total_loci,
		intrinsic_R,
		sigma_comp)

	genotypes_F = vcat([d.genotypes_F for d in pd.population]...)
	genotypes_M = vcat([d.genotypes_M for d in pd.population]...)
	x_locations_F = vcat([d.x_locations_F for d in pd.population]...)
	x_locations_M = vcat([d.x_locations_M for d in pd.population]...)
	y_locations_F = vcat([d.y_locations_F for d in pd.population]...)
	y_locations_M = vcat([d.y_locations_M for d in pd.population]...)



	@test count(i -> i == [0 0 0; 0 0 0], genotypes_F) == n^2 / 2
	@test count(i -> i == [0 0 0; 0 0 0], genotypes_M) == n^2 / 2
	@test count(i -> i == [1 1 1; 1 1 1], genotypes_F) == n^2 / 2
	@test count(i -> i == [1 1 1; 1 1 1], genotypes_M) == n^2 / 2
	@test length(genotypes_F) == n^2
	@test length(genotypes_M) == n^2
	@test length(x_locations_F) == n^2
	@test length(x_locations_M) == n^2
	@test length(y_locations_F) == n^2
	@test length(y_locations_M) == n^2

	[@test length(d.x_locations_F) == 1 for d in pd.population]
	[@test length(d.x_locations_M) == 1 for d in pd.population]
	[@test length(d.genotypes_F) == 1 for d in pd.population]
	[@test length(d.genotypes_M) == 1 for d in pd.population]
	[@test Population.assign_zone(pd.population[i].x_locations_F[1], pd.population[i].y_locations_F[1]) == i for i in eachindex(IndexCartesian(), pd.population)]
	[@test Population.assign_zone(pd.population[i].x_locations_M[1], pd.population[i].y_locations_M[1]) == i for i in eachindex(IndexCartesian(), pd.population)]


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

	for i in eachindex(genotypes)
		if x_locations[i] < 0.5
			@test genotypes[i] == [0 0 0; 0 0 0]
		else
			@test genotypes[i] == [1 1 1; 1 1 1]
		end
	end

	[@test length(pd.growth_rates_F[i]) == 1 for i in eachindex(IndexCartesian(), pd.population)]
end
@testset "new_location" begin
	locations = [Population.new_location([0.0f0, 1.0f0], [0.0f0, 1.0f0]) for i in 1:100000]
	function in_quadrant1(l)
		0 <= l[1] <= 0.5 && 0 <= l[2] <= 0.5
	end
	n = count(in_quadrant1, locations)

	@test abs(n / 100000 - 0.25) < 0.05
end

@testset "new_location_dispersal" begin
	x_locations = Float32[]
	y_locations = Float32[]
	sigma_disp = 0.01

	num_fail = 0

	for i in 1:1000
		x, y = Population.new_location(0.5f0, 0.5f0, sigma_disp)
		push!(x_locations, x)
		push!(y_locations, y)
		if x < 0.5 - 2 * sigma_disp || x > 0.5 + 2 * sigma_disp || y < 0.5 - 2 * sigma_disp || y > 0.5 + 2 * sigma_disp
			num_fail += 1
		end
	end

	average_x = sum(x_locations) / 1000

	average_y = sum(y_locations) / 1000

	@test num_fail < 70

	@test average_x > 0.495 && average_x < 0.505
	@test average_y > 0.495 && average_y < 0.505
end

@testset "genotype_mean" begin
	genotypes = [[1 0 1; 0 1 0], [0 0 1; 1 0 0], [1 1 1; 1 1 1]]

	@test Population.genotype_mean(genotypes[1], collect(1:2)) == 0.5

	@test Population.genotype_mean(genotypes[2], [2]) == 0
	@test Population.genotype_mean(genotypes[3], collect(1:3)) == 1
end

@testset "assign_zone" begin
	x_locations = [0.01f0, 0.32f0, 0.5f0]
	y_locations = [0.0f0, 0.75f0, 0.5f0]

	zones = Population.assign_zone.(x_locations, y_locations)

	@test zones == [CartesianIndex(1, 1), CartesianIndex(4, 8), CartesianIndex(6, 6)]
end

@testset "calc_survival_fitness_hetdisadvantage" begin
	genotypes = [[1 0 1; 0 1 0], [0 0 1; 0 0 1], [1 1 1; 1 1 1], [1 0 0; 0 1 0]]

	fitnesses = Population.calc_survival_fitness_hetdisadvantage.(genotypes, Ref(1))
	@test fitnesses == [1, 1, 1, 1]

	fitnesses = Population.calc_survival_fitness_hetdisadvantage.(genotypes, Ref(0))
	@test fitnesses == [0, 1, 1, 0]

	fitnesses = Population.calc_survival_fitness_hetdisadvantage.(genotypes, Ref(0.9))
	expected = [0.9, 1, 1, 0.932169]
	[@test abs(fitnesses[i] - expected[i]) < 0.001 for i in eachindex(expected)]
end

@testset "calc_survival_fitness_epistasis" begin
	genotypes = [[1 0 1; 0 1 0], [0 0 1; 0 0 1], [1 1 1; 1 1 1], [1 0 0; 0 1 0]]

	loci = 1:3
	fitnesses = Population.calc_survival_fitness_epistasis.(genotypes, Ref(loci), Ref(1))
	@test fitnesses == [1, 1, 1, 1]

	fitnesses = Population.calc_survival_fitness_epistasis.(genotypes, Ref(loci), Ref(0))
	@test fitnesses ≈ [0, 1 / 9, 1, 1 / 9]

	fitnesses = Population.calc_survival_fitness_epistasis.(genotypes, Ref(loci), Ref(0.9))
	expected = [0.9, 41 / 45, 1, 41 / 45]
	[@test abs(fitnesses[i] - expected[i]) < 0.001 for i in eachindex(expected)]
end

@testset "calc_species_overlap" begin
	pd1 = Population.PopulationData(10000, 2, 1.1, 0.01)
	population = pd1.population

	overlap, area_hybrid, area1, area2 = Population.calc_species_overlap(
		population,
		0.03,
		0.01,
		1:2,
	)

	@test overlap < 0.05
end
