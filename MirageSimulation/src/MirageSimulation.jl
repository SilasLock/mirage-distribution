module MirageSimulation
# This is the main simulation script. It's meant to be called from "init.jl" and nowhere else.
# Note that in order for Julia's package manager to run the project properly, this file needs
# to be named "MergeSimulation.jl". Otherwise, Julia claims the package/project isn't installed. 

# Grab solely the "Beta" and "cdf" functions from the Distributions package, and nothing else.
# Doing so avoids polluting the namespace.
using Distributions: Beta, cdf
# Do the same for the "dot" function from the LinearAlgebra package.
using LinearAlgebra: dot
# We're also going to use the Plots package; I think this should be all its relevant functions.
using Plots: plot, plot!, xlims!, title!, xlabel!, ylabel!
# TODO: Weirdly, "display()" isn't imported here but still seems to work? What's up with that?


# CDF section! You can specify new generators for
# the agent value distributions in here.
function kumaraswamyCDF(a::Float64, b::Float64)::Function
	@assert a > 0
	@assert b > 0
	return v -> max(min(1.0 - ((1.0 - (v^a))^(b - 1.0)), 1.0), 0.0)
end

function betaCDF(alpha::Float64, beta::Float64)::Function
    return v -> cdf(Beta(alpha, beta), v)
end

function mixtureZeroOneCDF(p::Float64)::Function
	return v -> {
		if v < 1.0
			return 1.0 - p
		else
			return 1.0
	}
end


# Dashboard section! You can code up new dashboards here and use them in a simulation.
function exampleDashboard(b::Float64)::Float64
	# This function is just random dashboard we're using as a a placeholder.
	# You'd code up an actual dashboard x(b) for most simulations instead of
	# this one.
	return max(min(1.0, 3 * b), 0.0)
end

function identityDashboard(b::Float64)::Float64
	# This dashboard just makes the probability that you get the item
	# equal to your bid. This is just about as simple as allocation
	# rules can get.
	@assert b >= 0
	@assert b <= 1
	return b
end

function straightlineDashboardFamily(b::Float64, theta::Float64)::Float64
	# The theta parameter is the ex ante allocation probability under a Unif([0, 1]) distribution of bids.
	@assert theta <= 1
	@assert theta >= 0
	@assert b >= 0
	@assert b <= 1
	if (theta <= 0.5)
		a = 2.0 * theta
		return a * b
	else
		a = (2.0 * theta) - 1.0
		return a + (1.0 - a) * b
	end
	return
end



function utilityWithIndices(v_index::Int64, b_index::Int64, x::Function, nonzerotypes::Int64)
	# This utility function should be correct and numerically stable for 64 bit floats.
	xarray = x.(collect(1:b_index) ./ nonzerotypes)
	return ((v_index - b_index) * x(b_index / nonzerotypes) + sum(xarray)) / nonzerotypes
end

function discretizedProbabilities(valueCDF::Function, nonzerotypes::Int64)
	# Using the CDF of the true value distribution F, calculate F(i / n) and
	# store the results in a giant array for all i in { 0, ..., n }.
	probabilities = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	# Then, starting from the top of the array and working down, subtract
	# F((i - 1) / n) from each value of F(i / n).
	pdfFromCDF!(probabilities)
	# for i = reverse(2:nonzerotypes + 1)
	# 	probabilities[i] -= probabilities[i - 1]
	# end
	return probabilities
end

function mirageCDF(v_hat::Float64, x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# The value of v_hat_index is "what is the largest type index i in {0,...,n} that
	# has a b_i greater than or equal to v_hat".
	v_hat_index::Int64 = floor(v_hat * nonzerotypes)
	# Construct the true pdf.
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	true_pdf_values = valueCDF.(agent_values)
	pdfFromCDF!(true_pdf_values)
	return mirageCDFWithIndices(v_hat_index, x, true_pdf_values, lambda, nonzerotypes)
end

function mirageCDFWithIndices(v_hat_index::Int64, x::Function, true_pdf_values::Vector{Float64}, lambda::Float64, nonzerotypes::Int64)::Float64
	# Calculates P(bid index <= `v_hat_index`).
	conditionalprobs = Vector{Float64}(undef, nonzerotypes + 1)
	for i in eachindex(conditionalprobs)
		# Compute the conditional cdf of the agent's bid, given that they have a true value index i - 1.
		conditionalprobs[i] = conditionalMirageCDFWithIndices(v_hat_index, i - 1, x, lambda, nonzerotypes)
	end
	# Calculate the probability of all true agent value regions.
	# probabilities = discretizedProbabilities(valueCDF, nonzerotypes)
	# Dotting the two vectors together calculates an expectation.
	return dot(conditionalprobs, true_pdf_values)
end

function mirageCDFImageWithIndices(x::Function, true_pdf_values::Vector{Float64}, lambda::Float64, nonzerotypes::Int64)::Vector{Float64}
	# Calculates P(bid index <= `v_hat_index`) for all possible values of `v_hat_index`.
	mirage_cdf_image = Vector{Float64}(undef, length(true_pdf_values))
	for j in eachindex(true_pdf_values)
		mirage_cdf_image[j] = mirageCDFWithIndices(j - 1, x, true_pdf_values, lambda, nonzerotypes)
	end
	return mirage_cdf_image
end

function conditionalMirageCDFWithIndices(v_hat_index::Int64, v_index::Int64, x::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# Calculates P(bid index <= `v_hat_index` | value index = `v_index`).
	if (v_hat_index >= v_index)
		return 1.0
	end
	weightarray = Vector{Float64}(undef, v_index + 1)
	for i = eachindex(weightarray)
		weightarray[i] = utilityWithIndices(v_index, i - 1, x, nonzerotypes)
	end
	weightarray = exp.(weightarray .* lambda)
	numerator = 0.0
	for i = 1:(v_hat_index + 1)
		numerator += weightarray[i]
	end
	denominator = sum(weightarray)
	return numerator / denominator
end

function conditionalMiragepdfWithIndices(v_hat_index::Int64, v_index::Int64, x::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# Calculates P(bid index = `v_hat_index` | value index = `v_index`).
	# Note that this should offer superior numerical stability compared
	# to taking consecutive differences of the `conditionalMirageCDFWithIndices` function.
	if (v_hat_index > v_index)
		return 0.0
	end
	weightarray = Vector{Float64}(undef, v_index + 1)
	for i = eachindex(weightarray)
		weightarray[i] = utilityWithIndices(v_index, i - 1, x, nonzerotypes)
	end
	weightarray = exp.(weightarray .* lambda)
	numerator = weightarray[v_hat_index + 1]
	denominator = sum(weightarray)
	return numerator / denominator
end

function exAnteAllocationProbability(x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# This is the probability q that dashboard x allocates an item to a quantal-responding agent.
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	true_pdf_values = valueCDF.(agent_values)
	pdfFromCDF!(true_pdf_values)
	mirage_pdf_values = mirageCDFImageWithIndices(x, true_pdf_values, lambda, nonzerotypes)
	pdfFromCDF!(mirage_pdf_values)
	# The following was the old way of generating the pdf from the CDF. Get rid of it at some point.
	# mirage_pdf_values = mirageCDFWithIndices.(collect(0:nonzerotypes), x, true_pdf_values, lambda, nonzerotypes)
	# for i = reverse(2:nonzerotypes + 1)
	# 	mirage_pdf_values[i] -= mirage_pdf_values[i - 1]
	# end
	dashboard_probs = x.(agent_values)
	return dot(dashboard_probs, mirage_pdf_values)
end

function pdfFromCDF!(probabilities::Vector{Float64})
	# Converts a vector representing a CDF into a vector representing the pdf.
	# Do this by taking the difference between adjacent CDF elements.
	for i = reverse(2:length(probabilities))
		probabilities[i] -= probabilities[i - 1]
	end
end


function inferValueDistribution(empiricalFrequency::Vector{Float64}, x::Function, lambda::Float64, nonzerotypes::Int64)::Vector{Float64}
	# TODO: Start here when you resume work next time! It's buggy for the zero type. =(
	matrixA = Matrix{Float64}(undef, nonzerotypes + 1, nonzerotypes + 1)
	for tempcartesianindex in CartesianIndices(matrixA)
		(j, i_hat) = Tuple(tempcartesianindex)
		matrixA[tempcartesianindex] = conditionalMiragepdfWithIndices(i_hat - 1, j - 1, x, lambda, nonzerotypes)
	end
	print("Made an A matrix: ")
	println(matrixA)
	onesVector = fill(1.0, nonzerotypes + 1)
	yVector = \(matrixA, onesVector)
	print("Made a yVector: ")
	println(yVector)
	matrixB = transpose(matrixA)
	for tempcartesianindex in CartesianIndices(matrixB)
		(i_hat, i) = Tuple(tempcartesianindex)
		matrixB[tempcartesianindex] *= yVector[i_hat]
	end
	inferredpdf = \(matrixB, empiricalFrequency)
	print("Made an inferred pdf: ")
	println(inferredpdf)
	return inferredpdf
end

function inferValueDistributionTest(empiricalFrequency::Vector{Float64}, x::Function, lambda::Float64, nonzerotypes::Int64)::Vector{Float64}
	@assert length(empiricalFrequency) == nonzerotypes + 1
	iterations = 100
	# TODO: This function is just a test of a particular iterative method! Let's see if it works properly before committing it.
	# fVector = Vector{Float64}(undef, nonzerotypes + 1)
	fVector = fill(1.0 / (nonzerotypes + 1), nonzerotypes + 1)
	muVector = fill(0.0, nonzerotypes + 1)
	gMatrix = Matrix{Float64}(undef, nonzerotypes + 1, nonzerotypes + 1)
	for tempcartesianindex in CartesianIndices(gMatrix)
		(i_hat, j) = Tuple(tempcartesianindex)
		gMatrix[tempcartesianindex] = conditionalMiragepdfWithIndices(i_hat - 1, j - 1, x, lambda, nonzerotypes)
	end
	intermediateMatrix = Matrix{Float64}(undef, nonzerotypes + 1, nonzerotypes + 1)
	gMatrixTranspose = transpose(gMatrix)
	# Each iteration, we
	# 1. map mu_i -> max(0, mu_i),
	# 2. solve for f from mu,
	# 3. map f_i -> max(0, f_i), then normalize sum f_i = 1,
	# 4. solve for mu from f.
	for t in 1:iterations
		# Step 1.
		for i in eachindex(muVector)
			muVector[i] = max(0.0, muVector[i])
			epsilon = 0.000000000001
			if (muVector[i] >= 1.0 - epsilon)
				muVector[i] = 1.0 - epsilon
			end
			if (fVector[i] > epsilon / (nonzerotypes + 1))
				muVector[i] *= muVector[i]
			end
		end
		# Step 2.
		tempVector = \(gMatrixTranspose, muVector)
		for tempcartesianindex in CartesianIndices(intermediateMatrix)
			(i_hat, j) = Tuple(tempcartesianindex)
			intermediateMatrix[tempcartesianindex] = gMatrix[tempcartesianindex] * (1.0 - tempVector[i_hat])
		end
		print("Mu vector: ")
		println(muVector)
		print("Temp vector: ")
		println(tempVector)
		print("Intermediate matrix: ")
		display(intermediateMatrix)
		println(intermediateMatrix[1, 2])
		print("Empirical frequency: ")
		println(empiricalFrequency)
		fVector = \(intermediateMatrix, empiricalFrequency)
		# Step 3.
		print("f vector before iteration: ")
		println(fVector)
		for i in eachindex(fVector)
			fVector[i] = max(0.0, fVector[i])
		end
		normalizer = 1.0 / sum(fVector)
		for i in eachindex(fVector)
			fVector[i] *= normalizer
		end
		print("f vector after iteration: ")
		println(fVector)
		# Step 4.
		tempVector = gMatrixTranspose * (empiricalFrequency ./ (gMatrix * fVector))
		for i in eachindex(muVector)
			muVector[i] = 1.0 - tempVector[i]
		end
		print("Mu vector after iteration: ")
		println(muVector)
	end
	return fVector
end


function uniformFrequency(nonzerotypes::Int64)::Vector{Float64}
	# Generates a uniform distribution's pdf, with zero weight on the zero type.
	# This function is useful for testing other functions that require an empirical frequency as input.
	uniformpdf = fill(1.0 / nonzerotypes, nonzerotypes + 1)
	# uniformpdf = Vector{Float64}(1.0 / nonzerotypes, nonzerotypes + 1)
	uniformpdf[1] = 0.0
	return uniformpdf
end

function generateSamples(numberofsamples::Int64, mirage_cdf_values::Vector{Float64})::Vector{Float64}
	# Generates `numberofsamples` from the mirage CDF, and returns them in the form of a frequency array.
	@assert numberofsamples >= 1
	empiricalCount = fill(0.0, length(mirage_cdf_values))
	for iter in 1:numberofsamples
		quantile = rand()
		for i in eachindex(mirage_cdf_values)
			# This implementation is extremely inefficient. Using binary search would be way faster.
			if (mirage_cdf_values[i] >= quantile)
				empiricalCount[i] += 1
				break
			end
		end
	end
	return empiricalCount ./ numberofsamples
end


function plotResults(values::Vector{Float64}, true_cdf_values::Vector{Float64}, mirage_cdf_values::Vector{Float64})
	x_axis_values = values
	y1_axis_values = true_cdf_values
	y2_axis_values = mirage_cdf_values

	ourplot = plot(x_axis_values, [y1_axis_values y2_axis_values], label=["True CDF" "Mirage CDF"], lw=[2 1])
	plot!(ourplot, legend=:outerbottom, legendcolumns=2)
	xlims!(ourplot, 0, 1)
	title!(ourplot, "True vs Mirage Distribution")
	xlabel!(ourplot, "v")
	ylabel!(ourplot, "F(v)")
	display(ourplot)
	println("Press ENTER when you're ready to stop looking at the plot.")
	junk = readline()
	# It would be really nice to exit a given graph after an arbitrary
	# keypress, but the tutorial here isn't as useful as one would like:
	# https://discourse.julialang.org/t/wait-for-a-keypress/20218/7
end

function inferExAnteAllocationProbability(x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# This is the probability q that dashboard x allocates an item to a quantal-responding agent.
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	true_pdf_values = valueCDF.(agent_values)
	pdfFromCDF!(true_pdf_values)
	mirage_pdf_values = mirageCDFImageWithIndices(x, true_pdf_values, lambda, nonzerotypes)
	pdfFromCDF!(mirage_pdf_values)
	# The following was the old way of generating the pdf from the CDF. Get rid of it at some point.
	# mirage_pdf_values = mirageCDFWithIndices.(collect(0:nonzerotypes), x, true_pdf_values, lambda, nonzerotypes)
	# for i = reverse(2:nonzerotypes + 1)
	# 	mirage_pdf_values[i] -= mirage_pdf_values[i - 1]
	# end
	dashboard_probs = x.(agent_values)
	return dot(dashboard_probs, mirage_pdf_values)
end

function testprintout(v::Float64, lambda::Float64, nonzerotypes::Int64)
	valueCDF = betaCDF(2.0, 2.0)
	print("\tValue v: ")
	println(v)
	print("True CDF at v = ")
	print(valueCDF(v))
	print("\tMirage CDF at v = ")
	println(mirageCDF(v, exampleDashboard, valueCDF, lambda, nonzerotypes))
end

function jokeyIntroSection()
	println("Simulation code will be run after this intro section.")
	println("For now, here's two outputs of our CDFs with a totally arbitrary parameterization I gave them.")
	print("Kumaraswamy distribution: ")
	println(kumaraswamyCDF(3.0, 6.0)(0.5))
	print("Beta distribution: ")
	println(betaCDF(2.0, 3.0)(0.6))
	print("Are these numbers probabilities between 0 and 1? Then great!")
	println(" We've already accomplished the first step of the project.")
	println("--------")
	println("To Aadityan: here are several printouts comparing the mirage and true CDF.")
	testprintout(0.7, 3.0, 50)
	testprintout(0.99, 3.0, 50)
	testprintout(0.2, 3.0, 50)
	println("Alright, let's do some simulations now!")
end

function displayMirage(x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)
	@assert nonzerotypes >= 1
	@assert lambda >= 0.0
	println("You've called the function to display both a true value CDF and its corresponding mirage CDF.")
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	true_cdf_values = valueCDF.(agent_values)
	true_pdf_values = valueCDF.(agent_values)
	pdfFromCDF!(true_pdf_values)
	mirage_cdf_values = mirageCDFImageWithIndices(x, true_pdf_values, lambda, nonzerotypes)
	print("Ex ante allocation probability (q): ")
	println(exAnteAllocationProbability(x, valueCDF, lambda, nonzerotypes))
	# If you were to instead use the command
	# mirage_cdf_values = mirageCDF.(agent_values, x, valueCDF, lambda, nonzerotypes)
	# it sometimes results in a floating point rounding error that causes some adjacent v_hat values to
	# produce the same output. And that's bad!
	plotResults(agent_values, true_cdf_values, mirage_cdf_values)
end

function startSimulation(numberofsamples::Int64, x_family::Function, theta_initial::Float64, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)
	@assert nonzerotypes >= 1
	@assert lambda >= 0.0
	println("Simulation Started!")
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	true_cdf_values = valueCDF.(agent_values)
	true_pdf_values = valueCDF.(agent_values)
	pdfFromCDF!(true_pdf_values)
	mirage_cdf_values = mirageCDFImageWithIndices(b -> x_family(b, theta_initial), true_pdf_values, lambda, nonzerotypes)
	empiricalFrequency = generateSamples(numberofsamples, mirage_cdf_values)
	print("Empirical frequency vector: ")
	println(empiricalFrequency)
	# empiricalFrequency = [2.0, 0.0, 0.1, 0.5, 0.2]
	inferred_value_pdf = inferValueDistributionTest(empiricalFrequency, b -> x_family(b, theta_initial), lambda, nonzerotypes)
	# Dot product inferred_value_pdf with the matrix of conditional cdfs of a *new* dashboard
	# return that

	theta_updated = theta_initial * theta_initial
	mirage_cdf_values_updated = mirageCDFImageWithIndices(b -> x_family(b, theta_updated), inferred_value_pdf, lambda, nonzerotypes)

	# print("Ex ante allocation probability (q): ")
	# println(exAnteAllocationProbability(b -> x_family(b, theta_initial), valueCDF, lambda, nonzerotypes))
	plotResults(agent_values, true_cdf_values, mirage_cdf_values_updated)
end


function main()
	jokeyIntroSection()
	x = identityDashboard
	valueCDF = betaCDF(2.0, 2.0)
	lambda = 30.0
	nonzerotypes = 90
	displayMirage(x, valueCDF, lambda, nonzerotypes)

	# Now let's test some inference procedures!
	# inferValueDistribution(uniformFrequency(nonzerotypes), x, lambda, nonzerotypes)
	nonzerotypes = 4
	# inferValueDistributionTest(uniformFrequency(nonzerotypes), x, lambda, nonzerotypes)
	numberofsamples = 10
	x_family = straightlineDashboardFamily
	theta_initial = 0.5
	startSimulation(numberofsamples, x_family, theta_initial, valueCDF, lambda, nonzerotypes)

	x = exampleDashboard
	valueCDF = betaCDF(2.0, 2.0)
	lambda = 30.0
	nonzerotypes = 90
	# displayMirage(x, valueCDF, lambda, nonzerotypes)
end

end # module MirageSimulation
