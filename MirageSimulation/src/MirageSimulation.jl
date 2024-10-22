module MirageSimulation
# This is the main simulation script. It's meant to be called from "init.jl" and nowhere else.
# Note that in order for Julia's package manager to run the project properly, this file needs
# to be named "MergeSimulation.jl". Otherwise, Julia claims the package/project isn't installed. 

# Grab solely the "Beta" and "cdf" functions from the Distributions package, and nothing else.
# Doing so avoids polluting the namespace.
using Distributions: Beta, cdf
# Do the same for the "dot" function from the LinearAlgebra package.
using LinearAlgebra: dot, mul!
# We're also going to use the Plots package; I think this should be all its relevant functions.
using Plots: plot, plot!, xlims!, title!, xlabel!, ylabel!, zlabel!, surface
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

# function mixtureZeroOneCDF(p::Float64)::Function
# 	return v -> {
# 		if v < 1.0
# 			return 1.0 - p
# 		else
# 			return 1.0
# 		end
# 	}
# end


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
		# Note that theta = 0.5 produces a diagonal line.
		a = (2.0 * theta) - 1.0
		return a + (1.0 - a) * b
	end
end

function postedpriceDashboardFamily(b::Float64, theta::Float64)::Float64
	# The theta parameter is the ex ante allocation probability under a Unif([0, 1]) distribution of bids.
	@assert theta <= 1
	@assert theta >= 0
	@assert b >= 0
	@assert b <= 1
	if (b > 1.0 - theta)
		return 1.0
	else
		return 0.0
	end
end

# function searchForThetaMatchingQ(x_family::Function, pdf_values::Vector{Float64}, q::Float64)
# 	@assert q <= 1
# 	@assert q >= 0
# 	# TODO: Implement this well! Remember that theta is a continuous parameter, you can make this work with that.
# 	# (Wait... actually, the fact that theta is continuous doesn't really help you here.)
# 	# When in doubt, err toward a higher theta. But also, don't err.
# 	# Maybe also print to the user how much they erred in their selection of theta due to discretization?
# 	return 0.5
# end


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
	if (isinf(lambda))
		if (v_hat_index == v_index)
			return 1.0
		else
			return 0.0
		end
	end
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

function conditionalMirageWeight(v_hat_index::Int64, v_index::Int64, x::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# Calculates the weight exp(lambda * utility).
	# This output is later used to calculate your mirage pdf.
	if (isinf(lambda))
		if (v_hat_index == v_index)
			return 1.0
		else
			return 0.0
		end
	end
	if (v_hat_index > v_index)
		return 0.0
	end
	# print("A test utility value: ")
	# println(utilityWithIndices(v_index, v_hat_index, x, nonzerotypes))
	return exp(lambda * utilityWithIndices(v_index, v_hat_index, x, nonzerotypes))
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

function inferValueDistributionFirstOrder(empiricalFrequency::Vector{Float64}, x::Function, lambda::Float64)::Vector{Float64}
	max_iterations = 1000
	supremumDistanceNeeded = 0.0000000000000001
	nonzerotypes = length(empiricalFrequency) - 1
	gVector = Vector{Float64}(undef, nonzerotypes + 1)
	gVectorFoundZero = Vector{Bool}(undef, nonzerotypes + 1)
	gMatrix = Matrix{Float64}(undef, nonzerotypes + 1, nonzerotypes + 1)
	intermediateMatrix = Matrix{Float64}(undef, nonzerotypes + 1, nonzerotypes + 1)
	for tempcartesianindex in CartesianIndices(gMatrix)
		(i_hat, j) = Tuple(tempcartesianindex)
		gMatrix[tempcartesianindex] = conditionalMiragepdfWithIndices(i_hat - 1, j - 1, x, lambda, nonzerotypes)
	end
	# We're going to acquire an initial estimate of fVector.
	fVector = \(gMatrix, empiricalFrequency)
	redistribute = 0.0
	for i in eachindex(fVector)
		# Sometimes, this initial estimate of fVector has negative components.
		# We're going to take those negative components and redistribute their total mass evenly to all components.
		if (fVector[i] < 0)
			redistribute -= fVector[i]
			fVector[i] = 0.0
		end
	end
	if (redistribute == 0)
		# Awesome! We got a valid probability vector on our first try.
		return fVector
	end
	# Looks like we didn't get the right vector the first time.
	# Let's iterate!
	supremumDistance = 0.0
	fVector .+= redistribute / length(fVector)
	piOverGVector = Vector{Float64}(undef, nonzerotypes + 1)
	for t in 1:max_iterations
		supremumDistance = 0.0
		gVectorFoundZero .= false
		mul!(gVector, gMatrix, fVector)
		for i_hat in eachindex(gVector)
			if (gVector[i_hat] == 0)
				gVectorFoundZero[i_hat] = true
				overwrite = 0.0
				for i in 1:length(empiricalFrequency)
					overwrite += gMatrix[i_hat, i]
				end
				gVector[i_hat] = overwrite
			end
		end
		piOverGVector .= empiricalFrequency ./ gVector
		for j in eachindex(fVector)
			new_f_j = 0.0
			for i in eachindex(piOverGVector)
				if (gVectorFoundZero[i])
					new_f_j += piOverGVector[i] * gMatrix[i, j]
				else
					new_f_j += piOverGVector[i] * gMatrix[i, j] * fVector[j]
				end
			end
			supremumDistance = max(supremumDistance, abs(fVector[j] - new_f_j))
			fVector[j] = new_f_j
		end
		# print("Supremum distance: ")
		# println(supremumDistance)
		if (supremumDistance < supremumDistanceNeeded)
			# println("Broke early!")
			break
		end
	end
	# println("Got a vector!")
	# print("Pi vector: ")
	# println(empiricalFrequency)
	# print("MLE f vector: ")
	# println(fVector)
	return fVector
end

# This is an old iterative method we're no longer using!
function inferValueDistributionTest(empiricalFrequency::Vector{Float64}, x::Function, lambda::Float64, nonzerotypes::Int64)::Vector{Float64}
	@assert length(empiricalFrequency) == nonzerotypes + 1
	# Change iterations to 100 to try to actually converge.
	# Use iterations=1 to do MLE without the "keep f(i) non-negative" constraint.
	iterations = 1
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


function plotResults(x_axis_label::AbstractString, y_axis_label::AbstractString, values::Vector{Float64}, first_cdf_name::AbstractString, first_cdf_values::Vector{Float64}, second_cdf_name::AbstractString, second_cdf_values::Vector{Float64})
	x_axis_values = values
	y1_axis_values = first_cdf_values
	y2_axis_values = second_cdf_values

	ourplot = plot(x_axis_values, [y1_axis_values y2_axis_values], label=[first_cdf_name second_cdf_name], lw=[2 1])
	plot!(ourplot, legend=:outerbottom, legendcolumns=2)
	xlims!(ourplot, 0, 1)
	title_name = "True vs Mirage Distribution"
	title_name = first_cdf_name * " vs " * second_cdf_name
	title!(ourplot, title_name)
	xlabel!(ourplot, x_axis_label)
	ylabel!(ourplot, y_axis_label)
	displayandpause(ourplot)
end

function visualizeErrors(theta_initials::Vector{Float64}, theta_updateds::Vector{Float64}, allocationproberrors::Array{Float64, 2})
	@assert length(theta_initials) == length(theta_updateds)
	# println(allocationproberrors)
	ourplot = surface(theta_initials, theta_updateds, allocationproberrors)
	# initialize a 3D plot with 1 empty series
	# ourplot = scatter3d(theta_initials, theta_updateds, allocationproberrors)
	# ourplot = plot3d(
	# 	1,
	# 	xlim = (-30, 30),
	# 	ylim = (-30, 30),
	# 	zlim = (0, 60),
	# 	title = "Lorenz Attractor",
	# 	legend = false,
	# 	marker = 2,
	# )
	title!(ourplot, "Errors in ex-ante allocation probability")
	xlabel!(ourplot, "theta old")
	ylabel!(ourplot, "theta new")
	zlabel!(ourplot, "Expected absolute error")
	displayandpause(ourplot)
end

function displayandpause(plotobject)
	display(plotobject)
	println("Press ENTER when you're ready to stop looking at the plot.")
	junk = readline()
	# It would be really nice to exit a given graph after an arbitrary
	# keypress, but the tutorial here isn't as useful as one would like:
	# https://discourse.julialang.org/t/wait-for-a-keypress/20218/7
end


function inferExAnteAllocationProbability(x::Function, mirage_CDF_values::Vector{Float64})::Float64
	# This is the probability q that dashboard x allocates an item to a quantal-responding agent.
	nonzerotypes = length(mirage_CDF_values) - 1
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	mirage_pdf_values = mirage_CDF_values
	pdfFromCDF!(mirage_pdf_values)
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
	plotResults("v", "F(v)", agent_values, "True CDF", true_cdf_values, "Mirage CDF", mirage_cdf_values)
end


function getExpectedQErrors(theta_granularity::Int64, numberofsamples::Int64, x_family::Function, true_pdf_values::Vector{Float64}, lambda::Float64, nonzerotypes::Int64)::Array{Float64, 2}
	@assert theta_granularity > 1
	# Observe that the error in q_hat, i.e. |q_overline_hat - q_overline|, is i.i.d sampled each iteration
	# It is also a random variable on the closed interval [0, 1]. Therefore, it has a maximum variance of 1/4.
	# This further implies that the variance of our sample average estimator of E[|q_overline_hat - q_overline|] is at most 1/4n
	# with n samples. Thus, the standard deviation is 0.5/sqrt(n).
	# If we pick iterations_to_compute_expected_error=10000, we therefore get a max standard deviation of 0.005, and that's great!
	iterations_to_compute_expected_error = 10000
	output_matrix = Array{Float64, 2}(undef, theta_granularity, theta_granularity)
	for j in 1:theta_granularity
		print("Progress: ")
		print(j - 1)
		print("/")
		println(theta_granularity)
		for i in 1:theta_granularity
			# Confirmed that in output_matrix[i, j], i refers to theta_updated and j refers to theta_initial.
			theta_initial = (j - 1) / (theta_granularity - 1)
			theta_updated = (i - 1) / (theta_granularity - 1)
			mirage_cdf_values = mirageCDFImageWithIndices(b -> x_family(b, theta_initial), true_pdf_values, lambda, nonzerotypes)
			mirage_cdf_values_updated_correct = mirageCDFImageWithIndices(b -> x_family(b, theta_updated), true_pdf_values, lambda, nonzerotypes)
			q_overline = inferExAnteAllocationProbability(b -> x_family(b, theta_updated), mirage_cdf_values_updated_correct)

			expected_error_metric = 0.0
			for expectation_iter in 1:iterations_to_compute_expected_error
				empiricalFrequency = generateSamples(numberofsamples, mirage_cdf_values)
				# print("Empirical frequency vector: ")
				# println(empiricalFrequency)
				inferred_value_pdf = inferValueDistributionFirstOrder(empiricalFrequency, b -> x_family(b, theta_initial), lambda)
				mirage_cdf_values_updated_inferred = mirageCDFImageWithIndices(b -> x_family(b, theta_updated), inferred_value_pdf, lambda, nonzerotypes)
				q_overline_hat = inferExAnteAllocationProbability(b -> x_family(b, theta_updated), mirage_cdf_values_updated_inferred)
				error_metric = abs(q_overline_hat - q_overline)
				expected_error_metric += error_metric
			end
			expected_error_metric /= iterations_to_compute_expected_error
			output_matrix[i, j] = expected_error_metric
		end
	end
	return output_matrix
end

function inferenceQErrorWithEpsEstimationGuarantee(theta_granularity::Int64, starting_theta::Float64, q_target::Float64, epsilon::Float64, x_family::Function, true_pdf_values::Vector{Float64}, lambda::Float64)
	println("Computing the (approximately) worst-case error in ex ante allocation probability you can get from using a new dashboard that intentionally targets q_target, when your inference from an old dashboard on the mirage CDF was at most epsilon in supremum distance.")
	nonzerotypes = length(true_pdf_values) - 1
	mirage_cdf_image = mirageCDFImageWithIndices(b -> x_family(b, starting_theta), true_pdf_values, lambda, nonzerotypes)
	mirage_cdf_image_upper = min.(mirage_cdf_image .+ epsilon, 1.0)
	mirage_cdf_image_lower = max.(mirage_cdf_image .- epsilon, 0.0)
	pdfFromCDF!(mirage_cdf_image_upper)
	pdfFromCDF!(mirage_cdf_image_lower)
	true_pdf_image_upper = inferValueDistributionFirstOrder(mirage_cdf_image_upper, b -> x_family(b, starting_theta), lambda)
	true_pdf_image_lower = inferValueDistributionFirstOrder(mirage_cdf_image_lower, b -> x_family(b, starting_theta), lambda)
	closest_theta_so_far_upper = 0.0
	smallest_q_error_so_far_upper = 100000000000.0
	closest_theta_so_far_lower = 0.0
	smallest_q_error_so_far_lower = 100000000000.0
	for i in 1:theta_granularity
		theta_check = (i - 1) / (theta_granularity - 1)
		mirage_cdf_image_upper_check = mirageCDFImageWithIndices(b -> x_family(b, theta_check), true_pdf_image_upper, lambda, nonzerotypes)
		q_upper_check = inferExAnteAllocationProbability(b -> x_family(b, theta_check), mirage_cdf_image_upper_check)
		if (abs(q_upper_check - q_target) < smallest_q_error_so_far_upper)
			closest_theta_so_far_upper = theta_check
			smallest_q_error_so_far_upper = abs(q_upper_check - q_target)
		end
		mirage_cdf_image_lower_check = mirageCDFImageWithIndices(b -> x_family(b, theta_check), true_pdf_image_lower, lambda, nonzerotypes)
		q_lower_check = inferExAnteAllocationProbability(b -> x_family(b, theta_check), mirage_cdf_image_lower_check)
		if (abs(q_lower_check - q_target) < smallest_q_error_so_far_lower)
			closest_theta_so_far_lower = theta_check
			smallest_q_error_so_far_lower = abs(q_lower_check - q_target)
		end
	end
	actual_mirage_cdf_image_upper = mirageCDFImageWithIndices(b -> x_family(b, closest_theta_so_far_upper), true_pdf_values, lambda, nonzerotypes)
	actual_ex_ante_upper = inferExAnteAllocationProbability(b -> x_family(b, closest_theta_so_far_upper), actual_mirage_cdf_image_upper)
	ex_ante_error_upper = abs(q_target - actual_ex_ante_upper)
	actual_mirage_cdf_image_lower = mirageCDFImageWithIndices(b -> x_family(b, closest_theta_so_far_lower), true_pdf_values, lambda, nonzerotypes)
	actual_ex_ante_lower = inferExAnteAllocationProbability(b -> x_family(b, closest_theta_so_far_lower), actual_mirage_cdf_image_lower)
	ex_ante_error_lower = abs(q_target - actual_ex_ante_lower)
	ex_ante_error = max(ex_ante_error_upper, ex_ante_error_lower)
	print("q target: ")
	println(q_target)
	print("epsilon (max supremum distance between estimated and true mirage CDF): ")
	println(epsilon)
	print("Ex ante error: ")
	println(ex_ante_error)
end


function startSimulation(theta_granularity::Int64, numberofsamples::Int64, x_family::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)
	@assert nonzerotypes >= 1
	@assert lambda >= 0.0
	iterations_to_compute_expected_error = 50

	println("Simulation Started!")
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	# true_cdf_values = valueCDF.(agent_values)
	true_pdf_values = valueCDF.(agent_values)
	pdfFromCDF!(true_pdf_values)

	x_thetas = collect(0:theta_granularity - 1) ./ (theta_granularity - 1)
	y_thetas = collect(0:theta_granularity - 1) ./ (theta_granularity - 1)

	z_errors = getExpectedQErrors(theta_granularity, numberofsamples, x_family, true_pdf_values, lambda, nonzerotypes)

	visualizeErrors(x_thetas, y_thetas, z_errors)
end

function calculateWelfare(theta_granularity::Int64, x_family::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)
	@assert nonzerotypes >= 1
	@assert lambda >= 0.0
	@assert theta_granularity > 1

	println("Starting to calculate welfare ratio.")
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	true_cdf_values = valueCDF.(agent_values)
	true_pdf_values = valueCDF.(agent_values)
	pdfFromCDF!(true_pdf_values)

	# x_thetas = collect(0:theta_granularity - 1) ./ (theta_granularity - 1)
	# y_thetas = collect(0:theta_granularity - 1) ./ (theta_granularity - 1)

	# This is a vector of the expectations of E[ x(B) | V] for all possible agent values V.
	conditional_allocation_probs = Vector{Float64}(undef, nonzerotypes + 1)

	# This is a vector of the probabilities P(B = b | V) for all possible bids b and a fixed V.
	conditional_mirage_pdf = Vector{Float64}(undef, nonzerotypes + 1)

	# This is a vector of x(B) for all possible agent bids B.
	xes = Vector{Float64}(undef, nonzerotypes + 1)

	# This is an adjustable vector detailing the optimal posted price mechanism's x() function.
	opt_xes = Vector{Float64}(undef, nonzerotypes + 1)

	# Arrays to store APX and OPT for various theta values in the APX dashboard.
	apx_at_theta = Vector{Float64}(undef, theta_granularity)
	opt_at_theta = Vector{Float64}(undef, theta_granularity)

	# What it says on the tin.
	ex_ante_allocation_prob_at_theta = Vector{Float64}(undef, theta_granularity)

	# Welfare ratio calculation.
	for t in 1:theta_granularity
		# print("Calculating ")
		# print(t)
		# print()
		theta = (t - 1) / (theta_granularity - 1)
		# Calculate the x(B) values ahead of time for this particular theta.
		for i_hat in eachindex(xes)
			xes[i_hat] = x_family(agent_values[i_hat], theta)
		end
		# Also reset the optimal posted price rule.
		for i in eachindex(opt_xes)
			opt_xes[i] = 1.0
		end
		# Calculate numerator.
		for j in eachindex(conditional_allocation_probs)
			# COMPLETED: This was rewritten to make it O(n^2) instead of O(n^3) in the number of discrete types!
			denominator = 0.0
			# denominator_array = fill(0.0, nonzerotypes + 1)
			for i_hat in 1:j
				# denominator_array[i_hat] = conditionalMirageWeight(i_hat - 1, j - 1, b -> x_family(b, theta), lambda, nonzerotypes)
				denominator += conditionalMirageWeight(i_hat - 1, j - 1, b -> x_family(b, theta), lambda, nonzerotypes)
			end
			# denominator = sum(denominator_array)
			for i_hat in eachindex(conditional_mirage_pdf)
				# conditional_mirage_pdf[i_hat] = conditionalMirageWeight(i_hat - 1, j - 1, b -> x_family(b, theta), lambda, nonzerotypes) / denominator
				conditional_mirage_pdf[i_hat] = conditionalMiragepdfWithIndices(i_hat - 1, j - 1, b -> x_family(b, theta), lambda, nonzerotypes)
			end
			# println(sum(conditional_mirage_pdf))
			conditional_allocation_probs[j] = dot(xes, conditional_mirage_pdf)
		end
		apx_at_theta[t] = sum(conditional_allocation_probs .* agent_values .* true_pdf_values)
		allocation_prob = dot(conditional_allocation_probs, true_pdf_values)
		ex_ante_allocation_prob_at_theta[t] = allocation_prob

		# Calculate denominator. This is the welfare from a (discretized) posted price mechanism with best response agents.
		# Also, the posted price mechanism has the same ex-ante allocation probability as the quantal response agents and their dashboard.
		# Note that we're doing this by hand, since the normal posted price dashboard family function doesn't know how to pick a theta
		# so they can allocate with probabilities between 0 and 1, which is necessary for this to work.

		# Find the first index i_star where F(i_star) >= 1 - q.
		i_star = nonzerotypes + 1
		for i in eachindex(true_cdf_values)
			opt_xes[i] = 0.0
			if (true_cdf_values[i] >= 1.0 - allocation_prob)
				i_star = i
				break
			end
		end
		# Do some clever math to find the denominator.
		denominator = sum(agent_values .* opt_xes .* true_pdf_values) + (agent_values[i_star] * (true_cdf_values[i_star] - 1.0 + allocation_prob))
		opt_at_theta[t] = denominator
	end
	# println(apx_at_theta)
	# println(opt_at_theta)
	welfare_ratios = Vector{Float64}(undef, theta_granularity)
	for i in eachindex(welfare_ratios)
		if (opt_at_theta[i] == 0.0)
			# If the optimal welfare is 0, then we know our apx welfare (which can't do better) must also be zero.
			# In these cases, we say that we've captured 100% of the welfare, and get an approximation ratio of 1.
			welfare_ratios[i] = 1.0
		else
			# Otherwise, compute the approximation ratio the normal way.
			welfare_ratios[i] = apx_at_theta[i] / opt_at_theta[i]
		end
	end
	plotResults("theta", "APX/OPT", collect(0:theta_granularity - 1) ./ (theta_granularity - 1), "Welfare ratio", welfare_ratios, "Ex ante allocation probability", ex_ante_allocation_prob_at_theta)
end


function main()
	jokeyIntroSection()
	x = identityDashboard
	# valueCDF = betaCDF(2.0, 2.0)
	# valueCDF = betaCDF(100.0, 3.0)
	valueCDF = betaCDF(2.0, 1000.0)
	lambda = 30.0
	nonzerotypes = 90
	# displayMirage(x, valueCDF, lambda, nonzerotypes)


	# Now let's test some inference procedures!
	# inferValueDistribution(uniformFrequency(nonzerotypes), x, lambda, nonzerotypes)
	nonzerotypes = 10
	# inferValueDistributionTest(uniformFrequency(nonzerotypes), x, lambda, nonzerotypes)
	theta_granularity = 100
	numberofsamples = 100
	x_family = straightlineDashboardFamily
	# x_family = postedpriceDashboardFamilys
	# startSimulation(theta_granularity, numberofsamples, x_family, valueCDF, lambda, nonzerotypes)
	# calculateWelfare(theta_granularity, x_family, valueCDF, lambda, nonzerotypes)

	theta_granularity = 100
	starting_theta = 0.5
	q_target = 0.3
	epsilon = 0.02
	x_family = straightlineDashboardFamily
	valueCDF = betaCDF(2.0, 10.0)
	true_pdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	pdfFromCDF!(true_pdf_values)
	lambda = 3.0
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	epsilon = 0.2
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)



	# x = exampleDashboard
	# valueCDF = betaCDF(2.0, 2.0)
	# lambda = 30.0
	# nonzerotypes = 90
	# displayMirage(x, valueCDF, lambda, nonzerotypes)
end

end # module MirageSimulation
