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
using Plots: plot, plot!, xlims!, title!, xlabel!, ylabel!, zlabel!, surface, savefig
# TODO: Weirdly, "display()" isn't imported here but still seems to work? What's up with that?
# Acquire the ability to do weighted categorical sampling from StatsBase.
using StatsBase: sample, Weights
# Try to import the Julia profiler for use with the @profile macro.
using Profile

# This is a temporary hack to stop pausing on image displays when using a remote machine.
not_using_a_remote_machine = false

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

function pointMassCDF(pointmass::Float64)::Function
	return v ->
		if (v < pointmass)
			return 0.0
		else
			return 1.0
		end
end

function twoPointMassesCDF(pointmass_low::Float64, pointmass_high::Float64, prob_high::Float64)::Function
	@assert pointmass_high >= pointmass_low
	return v ->
	if (v < pointmass_low)
		return 0.0
	elseif (v < pointmass_high)
		return 1.0 - prob_high
	else
		return 1.0
	end
end

function mixtureZeroOneCDF(p::Float64, lowest::Float64)::Function
	# p is the probability of getting value 1, lowest is a number \approx 0 which occurs with probability 1 - p.
	return v -> 
		if v < lowest
			return 0.0
		elseif v < 0.9999
			return 1.0 - p
		else
			return 1.0
		end
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

function splineSigmoidGenerator(p_high::Float64)::Function
	f = (b::Float64, theta::Float64) -> begin
		# The theta parameter is 1, minus the ratio of p_low to p_high.
		@assert theta <= 1
		@assert theta >= 0
		@assert b >= 0
		@assert b <= 1
		p_low = (1.0 - theta) * p_high
		normed_b = (b - p_low) / (p_high - p_low)
		if (b < p_low)
			return 0.0
		elseif (b < (p_low + p_high) / 2.0)
			return 2 * (normed_b * normed_b)
		elseif (b < p_high)
			return -1.0 + 2.0 * normed_b * (2.0 - normed_b)
		else
			return 1.0
		end
	end
	return f
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

# TODO: Left off here! This new pdf image calculator will be much more performant than the CDF image calculator.
# It's complete, but you should update the rest of the codebase to use it, since it's so much more performant than `mirageCDFImageWithIndices`.
function miragePDFImageWithIndices(x::Function, true_pdf_values::Vector{Float64}, lambda::Float64)::Vector{Float64}
	# Calculates P(bid index == `v_hat_index`) for all possible values of `v_hat_index`.
	if (isinf(lambda))
		# Assume that we're using a truthful dashboard.
		return copy(true_pdf_values)
	end
	nonzerotypes = length(true_pdf_values) - 1
	mirage_pdf_image = Vector{Float64}(undef, length(true_pdf_values))
	reweighted_value_pdf = Vector{Float64}(undef, length(true_pdf_values))
	for j in eachindex(reweighted_value_pdf)
		denominator::Float64 = 0.0
		for i = 1:j
			denominator += exp(lambda * utilityWithIndices(j - 1, i - 1, x, nonzerotypes))
		end
		reweighted_value_pdf[j] = true_pdf_values[j] / denominator
	end
	weightarray::Vector{Float64} = fill(0.0, length(true_pdf_values))
	for j in eachindex(mirage_pdf_image)
		for i in j:length(mirage_pdf_image)
			weightarray[i] = exp(lambda * utilityWithIndices(i - 1, j - 1, x, nonzerotypes))
		end
		# println(weightarray)
		# println()
		mirage_pdf_image[j] = dot(weightarray, reweighted_value_pdf)
		weightarray[j] = 0.0
	end
	return mirage_pdf_image
end

function conditionalMirageCDFWithIndices(v_hat_index::Int64, v_index::Int64, x::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# Calculates P(bid index <= `v_hat_index` | value index = `v_index`).
	if (v_hat_index >= v_index)
		return 1.0
	elseif (isinf(lambda))
		return 0.0
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
	# This algorithm is basically a variant of the EM algorithm, and it has quite a poor convergence rate.
	max_iterations = 10000
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
	if (supremumDistance >= supremumDistanceNeeded)
		println("WARNING: Iterative method failed to converge.")
		print("PDF was still moving by size ")
		print(supremumDistance)
		print(" when we needed below ")
		println(supremumDistanceNeeded)
	end
	# println("Got a vector!")
	# print("Pi vector: ")
	# println(empiricalFrequency)
	# print("MLE f vector: ")
	# println(fVector)
	return fVector
end

# This is an old iterative method we're no longer using!
function inferValueDistributionTest(empiricalFrequency::Vector{Float64}, x::Function, lambda::Float64)::Vector{Float64}
	# @assert length(empiricalFrequency) == nonzerotypes + 1
	nonzerotypes = length(empiricalFrequency) - 1
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
	if (not_using_a_remote_machine)
		display(plotobject)
		println("Press ENTER when you're ready to stop looking at the plot.")
		junk = readline()
	end
	# It would be really nice to exit a given graph after an arbitrary
	# keypress, but the tutorial here isn't as useful as one would like:
	# https://discourse.julialang.org/t/wait-for-a-keypress/20218/7
end


function inferExAnteAllocationProbabilityFromCDF(x::Function, mirage_CDF_values::Vector{Float64})::Float64
	# This is the probability q that dashboard x allocates an item to a quantal-responding agent.
	nonzerotypes = length(mirage_CDF_values) - 1
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	mirage_pdf_values = copy(mirage_CDF_values)
	pdfFromCDF!(mirage_pdf_values)
	dashboard_probs = x.(agent_values)
	return dot(dashboard_probs, mirage_pdf_values)
end

function inferExAnteAllocationProbability(x::Function, mirage_PDF_values::Vector{Float64})::Float64
	# This is the probability q that dashboard x allocates an item to a quantal-responding agent.
	nonzerotypes = length(mirage_PDF_values) - 1
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	dashboard_probs = x.(agent_values)
	return dot(dashboard_probs, mirage_PDF_values)
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
			q_overline = inferExAnteAllocationProbabilityFromCDF(b -> x_family(b, theta_updated), mirage_cdf_values_updated_correct)

			expected_error_metric = 0.0
			for expectation_iter in 1:iterations_to_compute_expected_error
				empiricalFrequency = generateSamples(numberofsamples, mirage_cdf_values)
				# print("Empirical frequency vector: ")
				# println(empiricalFrequency)
				inferred_value_pdf = inferValueDistributionFirstOrder(empiricalFrequency, b -> x_family(b, theta_initial), lambda)
				mirage_cdf_values_updated_inferred = mirageCDFImageWithIndices(b -> x_family(b, theta_updated), inferred_value_pdf, lambda, nonzerotypes)
				q_overline_hat = inferExAnteAllocationProbabilityFromCDF(b -> x_family(b, theta_updated), mirage_cdf_values_updated_inferred)
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
	mirage_cdf_image_upper[length(mirage_cdf_image_upper)] = 1.0
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
		q_upper_check = inferExAnteAllocationProbabilityFromCDF(b -> x_family(b, theta_check), mirage_cdf_image_upper_check)
		if (abs(q_upper_check - q_target) < smallest_q_error_so_far_upper)
			closest_theta_so_far_upper = theta_check
			smallest_q_error_so_far_upper = abs(q_upper_check - q_target)
		end
		mirage_cdf_image_lower_check = mirageCDFImageWithIndices(b -> x_family(b, theta_check), true_pdf_image_lower, lambda, nonzerotypes)
		q_lower_check = inferExAnteAllocationProbabilityFromCDF(b -> x_family(b, theta_check), mirage_cdf_image_lower_check)
		if (abs(q_lower_check - q_target) < smallest_q_error_so_far_lower)
			closest_theta_so_far_lower = theta_check
			smallest_q_error_so_far_lower = abs(q_lower_check - q_target)
		end
	end
	actual_mirage_cdf_image_upper = mirageCDFImageWithIndices(b -> x_family(b, closest_theta_so_far_upper), true_pdf_values, lambda, nonzerotypes)
	actual_ex_ante_upper = inferExAnteAllocationProbabilityFromCDF(b -> x_family(b, closest_theta_so_far_upper), actual_mirage_cdf_image_upper)
	ex_ante_error_upper = abs(q_target - actual_ex_ante_upper)
	actual_mirage_cdf_image_lower = mirageCDFImageWithIndices(b -> x_family(b, closest_theta_so_far_lower), true_pdf_values, lambda, nonzerotypes)
	actual_ex_ante_lower = inferExAnteAllocationProbabilityFromCDF(b -> x_family(b, closest_theta_so_far_lower), actual_mirage_cdf_image_lower)
	ex_ante_error_lower = abs(q_target - actual_ex_ante_lower)
	if (ex_ante_error_upper > ex_ante_error_lower)
		actual_theta_new = closest_theta_so_far_upper
	else
		actual_theta_new = closest_theta_so_far_lower
	end
	ex_ante_error = max(ex_ante_error_upper, ex_ante_error_lower)
	print("Initial theta value: ")
	println(starting_theta)
	print("New theta value to achieve q target: ")
	println(actual_theta_new)
	print("q target: ")
	println(q_target)
	print("epsilon (max supremum distance between estimated and true mirage CDF): ")
	println(epsilon)
	print("Ex ante error: ")
	println(ex_ante_error)
	# print("Ex ante error divided by q target")
	println("------------------")
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

function calculateWelfareOfFamily(theta_granularity::Int64, x_family::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)
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


function testingTheSigmoidGenerator()
	lambda = 0.15

	p_high = 0.3
	chosenSpline = splineSigmoidGenerator(p_high)
	println("Here's the dashboard:")
	for i = 0:100
		println(chosenSpline(1.0 - (i / 100.0), 0.2))
	end
	println(chosenSpline(0.4, 0.2))
	println("^^ This is the dashboard.")
	# TODO: This is where we left off!
	println("Press Enter to continue.")
	junk = readline()
end

# function getThetaAssociatedWithQ(target_q::Float64, theta_granularity::Int64, dashboard_family::Function, value_pdf::Vector{Float64}, lambda::Float64)::Float64
# 	# TODO: This seems to be a duplicate function of findThetaWhichObtainsQ! You might want to deduplicate this code.
# 	# This version seems better than the other one, because of its usage of theta_granularity as an explicit parameter.
# 	# But it also doesn't output a tight lower bound, which might be less good?
#	# TODO: Delete this function! It's old, glitchy, and bad. It's entirely superseded by findThetaWhichObtainsQ.

# 	# IMPORTANT: This function assumes that the allocation probability is *increasing* in theta!
# 	# Dashboard families must be defined so that this assumption holds true.
# 	theta_index_high = theta_granularity
# 	theta_index_low = 0

# 	while (theta_index_high > theta_index_low)
# 		theta = (theta_index_high + theta_index_low) / (2.0 * theta_granularity)
# 		mirage_cdf = mirageCDFImageWithIndices(b -> dashboard_family(b, theta), value_pdf, lambda, length(value_pdf) - 1)
# 		mid_q = inferExAnteAllocationProbabilityFromCDF(b -> dashboard_family(b, theta), mirage_cdf)
# 		if (mid_q > target_q)
# 			theta_index_high = round((theta_index_high + theta_index_low) / 2.0)
# 		else
# 			theta_index_low = round((theta_index_high + theta_index_low) / 2.0)
# 		end
# 	end
# 	# mirage_cdf_high = mirageCDFImageWithIndices(b -> dashboard_family(b, theta_high), value_pdf, lambda, length(value_pdf) - 1)
# 	# mirage_cdf_low = mirageCDFImageWithIndices(b -> dashboard_family(b, theta_low), value_pdf, lambda, length(value_pdf) - 1)

# 	return theta
# end

function getIndexSamplesFromPdf(pdf_to_be_sampled, num_samples::Int64)
	# Sample num_samples indices from a distribution. (Usually the true value distribution or the mirage distribution.)
	# TODO: Make this take an rng parameter, like Xoshiro!
	return sample(1:length(pdf_to_be_sampled), Weights(pdf_to_be_sampled), num_samples)
end

function frequencyFromIndexSamples(sample_indices::Array{Int64}, nonzerotypes)::Array{Float64}
	# Convert a index samples vector into an empirical frequency vector (i.e. a PDF).
	countVector = fill(0.0, nonzerotypes + 1)
	for i in sample_indices
		countVector[i] += 1
	end
	# println(countVector ./ length(sample_indices))
	return countVector ./ length(sample_indices)
	# return fill(1.0 / (nonzerotypes + 1), nonzerotypes + 1)
end

function findThetaWhichObtainsQ(x_family::Function, q_target::Float64, value_pdf::Vector{Float64}, lambda::Float64)::Float64
	# Uses binary search to find a parameter theta of a dashboard family which gives allocation probability q_target.
	# If there doesn't exist an exact theta value which provides this exact allocation probability, we output something *right* below it.
	# That way, ex-ante supply constraints are always satisfied, even if they aren't perfectly tight.

	nonzerotypes = length(value_pdf) - 1
	# Hard code the theta granularity, which determines how carefully we search through parameter space, to be
	# at least the number of non-zero agent types times the maximal constant multiplier we can use.
	# We compute the highest multiplier with integer division.
	highest_integer_multiplier = div(typemax(Int64), nonzerotypes)
	# ... then we multiply it by the number of non-zero agent types.
	theta_granularity = highest_integer_multiplier * nonzerotypes
	# This gives us an integer number of theta values which can fall in between different thetas equal to the type.
	# This is particularly elegant for computing thetas from the posted price family (assuming that the tie-breaking rules work right.)
	@assert theta_granularity > 0

	# IMPORTANT: This function assumes that the allocation probability is *increasing* in theta!
	# Dashboard families must be defined so that this assumption holds true.

	# theta_index_low and theta_index_high are the boundaries of the binary search region. They get closer to one another as we proceed.
	theta_index_low = 0
	theta_index_high = theta_granularity

	while (theta_index_low < theta_index_high)
		# theta_index_mid is where we target for testing. It's designed to both (a) round up and (b) to avoid integer overflow.
		# println("Divs calculating...")
		theta_index_mid = theta_index_high - div(theta_index_high - theta_index_low, 2)
		theta_guess::Float64 = theta_index_mid / theta_granularity
		# println("Curried function calculating...")
		x = b -> x_family(b, theta_guess)


		# print("Theta low: ")
		# println(theta_index_low / theta_granularity)
		# print("Theta high: ")
		# println(theta_index_high / theta_granularity)

		# println("Mirage CDF calculating...")
		# mirage_cdf = mirageCDFImageWithIndices(x, value_pdf, lambda, nonzerotypes)
		# println(mirage_cdf)
		# println("q_guess calculating...")
		# q_guess = inferExAnteAllocationProbabilityFromCDF(x, mirage_cdf)

		# println("Mirage PDF calculating...")
		mirage_pdf = miragePDFImageWithIndices(x, value_pdf, lambda)
		# println("q_guess calculating...")
		q_guess = inferExAnteAllocationProbability(x, mirage_pdf)

		# mirage_cdf_to_pdf = copy(mirage_cdf)
		# pdfFromCDF!(mirage_cdf_to_pdf)
		# println(mirage_cdf_to_pdf)
		# println(maximum(abs.(mirage_cdf_to_pdf - mirage_pdf)))
		# println(sum(mirage_cdf_to_pdf))
		# println(sum(mirage_pdf))
		# println(sum(value_pdf))
		

		# print("Q guess: ")
		# println(q_guess)
		# print("Q target: ")
		# println(q_target)

		# Some notes: This function will successfully return the highest theta that produces a q_guess which satisfies the ex-ante supply constraint.
		# This even works in cases where there is no precise theta which gives us the desired ex-ante allocation probability,
		# or in cases where there are multiple theta that produce allocation probabilities exactly equal to q_target.
		# In particular, it will always output theta = 1.0 when q_target was set to 1.0. However, it may not always output theta = 0.0 when q_target = 0.0.
		# If there are no sufficiently low parameters theta which allow you to satisfy a low q_target, the function returns theta = 0.0.
		if (q_guess <= q_target)
			theta_index_low = theta_index_mid
		else
			theta_index_high = theta_index_mid - 1
		end
	end
	# println(theta_index_low / theta_granularity)

	# Return a (hopefully very tight) lower bound on the correct theta parameter.
	return theta_index_low / theta_granularity
end

function calculateWelfare(x::Function, true_pdf_values::Vector{Float64}, lambda::Float64)::Float64
	welfare_counter = 0.0
	nonzerotypes = length(true_pdf_values) - 1
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	for i in eachindex(agent_values)
		conditional_allocation = 0.0
		for j in eachindex(agent_values)
			conditional_allocation += x(agent_values[j]) * conditionalMiragepdfWithIndices(j - 1, i - 1, x, lambda, nonzerotypes)
		end
		welfare_counter += true_pdf_values[i] * agent_values[i] * conditional_allocation
	end
	return welfare_counter
end

function sigmoidEstimationProcedure(true_cdf_values::Vector{Float64}, q::Float64, w::Float64, lambda::Float64, num_samples::Int64, num_rounds::Int64)
	# TODO: Make these better initial estimated parameters! Right now they're pretty randomly chosen. You do want them hard-coded, just...
	# ... not quite these initial estimates.
	p_high_initial = 0.7
	p_low_initial = 0.2

	true_pdf_values = copy(true_cdf_values)
	pdfFromCDF!(true_pdf_values)

	nonzerotypes = length(true_pdf_values) - 1

	# Note: Normally, for the 1/2 * 3/4 approximation, w should equal 1/2.

	# Keep track of inference error.
	sum_of_q_abs_errors = 0.0
	# Keep track of the fraction of the time we select an ex-ante allocation probability that is too high.
	num_rounds_too_high_exante_allocation_prob = 0
	# Keep track of the total overallocation probability across all rounds.
	cumulative_overallocation_probability = 0.0

	# Keep track of sums of apx welfare and opt welfare.
	# We plan to sum these across all rounds, and then divide one by the other.
	sum_of_apx_welfares = 0.0
	# optimal_br_postedprice = b -> postedpriceDashboardFamily(b, findThetaWhichObtainsQ(postedpriceDashboardFamily, q, true_pdf_values, Inf))
	opt_br_welfare = 0.0
	# sanity_check = 0.0
	for i in reverse(eachindex(true_cdf_values))
		# println(sanity_check)
		F_iMinus1::Float64 = (i > 1) ? true_cdf_values[i - 1] : 0.0
		if (q < 1.0 - F_iMinus1)
			# q >= 1.0 - true_cdf_values[i]
			# q == 1.0 - true_cdf_values[i] + tiny_alloc_prob
			# tiny_alloc_prob = q - 1.0 + true_cdf_values[i]
			tiny_alloc_prob = q - 1.0 + true_cdf_values[i]
			opt_br_welfare += (i - 1) * tiny_alloc_prob
			# sanity_check += tiny_alloc_prob
			break
		end
		opt_br_welfare += (i - 1) * true_pdf_values[i]
		# sanity_check += true_pdf_values[i]
	end
	# print("Sanity check: ")
	# println(sanity_check)
	opt_br_welfare /= nonzerotypes
	# opt_br_welfare = calculateWelfare(optimal_br_postedprice, true_pdf_values, Inf)
	print("Optimal welfare: ")
	println(opt_br_welfare)

	print("Target q: ")
	println(q)

	print("w: ")
	println(w)

	print("lambda: ")
	println(lambda)


	previous_spline = b -> splineSigmoidGenerator(p_high_initial)(b, 1.0 - (p_low_initial / p_high_initial))

	for i in 1:num_rounds
		# Use the sigmoid from the last round (t - 1) to get n samples that we use to estimate the value distribution.
		mirage_pdf_image = miragePDFImageWithIndices(previous_spline, true_pdf_values, lambda)
		# mirage_pdf_image = mirageCDFImageWithIndices(previous_spline, true_pdf_values, lambda, nonzerotypes)
		# pdfFromCDF!(mirage_pdf_image)
		empirical_samples = getIndexSamplesFromPdf(mirage_pdf_image, num_samples)
		# TODO: Confirm that the following line works and is not buggy!
		empirical_frequency_vector = frequencyFromIndexSamples(empirical_samples, nonzerotypes)
		# Use the samples from the mirage to infer the maximum-likelihood estimated value pdf.
		# TODO: Confirm that this estimation method is functional!
		estimated_value_pdf = inferValueDistributionFirstOrder(empirical_frequency_vector, previous_spline, lambda)
		# estimated_value_pdf = copy(mirage_pdf_image)

		# We now have an estimate of the value distribution.
		# Assuming that this estimate is correct, proceed to the steps in the following section:

		
		# Find out what p_low and p_high are needed to achieve a particular allocation probability of a posted price dashboard.
		# Use a posted price dashboard, do binary search to find the price p_high
		# at which the allocation prob is wq where w is a constant in [0, 1]
		# Note that we use 1 - theta to calculate p_high, since the theta parameter for a posted price dashboard is "1 minus the price."
		p_high = 1.0 - findThetaWhichObtainsQ(postedpriceDashboardFamily, w * q, estimated_value_pdf, lambda)
		print("P high: ")
		println(p_high)
		# Note: We have a proof that binary search should be successful for locating the right theta to obtain allocation probability q,
		# which works by showing that the allocation probability is monotone increasing in theta.

		# Use a sigmoid dashboard, which we make the sigmoid dashboard a curried function of p_high, to find the price p_low which
		# makes the overall sigmoid dashboard have ex-ante allocation probability q.
		lower_theta = findThetaWhichObtainsQ(splineSigmoidGenerator(p_high), q, estimated_value_pdf, lambda)
		p_low = p_high * (1.0 - lower_theta)
		print("P low: ")
		println(p_low)

		# Note that both these simulations to find p_high and p_low are not run on actual agents. Instead, they're run on the inferred value distribution we
		# constructed from the data collected by a sigmoid dashboard in the previous timestep (t - 1).

		# Construct the sigmoid dashboard from p_high and p_low (i.e. lower_theta).
		chosenSpline = b -> splineSigmoidGenerator(p_high)(b, lower_theta)

		# What do we do with a sigmoid dashboard now?
		# Sanity check: make sure the sigmoid with p_low and p_high has an allocation probability in between
		# q and wq.
		mirage_pdf_image = miragePDFImageWithIndices(chosenSpline, true_pdf_values, lambda)
		allocation_prob = inferExAnteAllocationProbability(chosenSpline, mirage_pdf_image)
		# mirage_cdf_image = mirageCDFImageWithIndices(chosenSpline, true_pdf_values, lambda, nonzerotypes)
		# allocation_prob = inferExAnteAllocationProbabilityFromCDF(chosenSpline, mirage_cdf_image)
		print("Allocation probability: ")
		println(allocation_prob)

		# If we have an ex-ante allocation probability higher than our supply constraint/target, do an adjustment where we
		# cancel the auction entirely with some probability, getting 0 welfare when this occurs.
		# This is a brute-force, sledgehammer-ey way to force the ex-ante allocation probability below our target.
		adjustment = (q < allocation_prob) ? (q / allocation_prob) : 1.0

		# Calculate/tally sum of welfares of spline sigmoid dashboard.
		current_welfare = adjustment * calculateWelfare(chosenSpline, true_pdf_values, lambda)
		print("Current welfare: ")
		println(current_welfare)
		print("Optimal welfare: ")
		println(opt_br_welfare)
		print("Welfare ratio compared to BR at allocation probability q: ")
		println(current_welfare / opt_br_welfare)
		sum_of_apx_welfares += current_welfare

		# Calculate/tally inference error.
		sum_of_q_abs_errors += abs(allocation_prob - q)
		print("Absolute inference error: ")
		println(abs(allocation_prob - q))

		if (allocation_prob > q)
			print("A bad event occurred! We overallocated by a probability margin of ")
			println(allocation_prob - q)
		end

		# Update the counter of how many times the ex-ante allocation probability has been above our target.
		num_rounds_too_high_exante_allocation_prob += (allocation_prob > q)

		# If we had a too-high ex-ante allocation probability, record the fraction of agents that would have been allocated but now can't.
		# We're keeping track of the sum of these values to compute their average value across rounds at the end.
		cumulative_overallocation_probability += max(0.0, allocation_prob - q)
		# cumulative_overallocation_probability += (allocation_prob > q) * (allocation_prob - q)

		# Use sigmoid dashboard in next round (t) to get another n samples, and get another estimate of the value distribution.
		previous_spline = chosenSpline
	end

	sum_of_opt_br_welfares = num_rounds * opt_br_welfare
	# print("Averaged welfare ratio: ")
	# println(sum_of_apx_welfares / sum_of_opt_br_welfares)

	# print("Worst-case welfare ratio assuming that we fit p_high correctly: ")
	# println(w * 3/4)

	# print("Average absolute inference error: ")
	# println(sum_of_q_abs_errors / num_rounds)

	# print("Fraction of rounds where the ex-ante allocation probability was too high: ")
	# println(num_rounds_too_high_exante_allocation_prob / num_rounds)

	# print("Expected fraction of agents who were meant to be allocated among the q fraction, but couldn't be because we made the ex-ante allocation probability too high: ")
	# println(cumulative_overallocation_probability / num_rounds)
	
	# mirage_cdf_image_upper = min.(mirage_cdf_image .+ epsilon, 1.0)
	# mirage_cdf_image_upper[length(mirage_cdf_image_upper)] = 1.0
	# mirage_cdf_image_lower = max.(mirage_cdf_image .- epsilon, 0.0)
	# pdfFromCDF!(mirage_cdf_image_upper)
	# pdfFromCDF!(mirage_cdf_image_lower)
	# true_pdf_image_upper = inferValueDistributionFirstOrder(mirage_cdf_image_upper, b -> x_family(b, starting_theta), lambda)
	return SigmoidSummary(
		sum_of_apx_welfares / sum_of_opt_br_welfares,
		w * 3/4,
		sum_of_q_abs_errors / num_rounds,
		num_rounds_too_high_exante_allocation_prob / num_rounds,
		cumulative_overallocation_probability / num_rounds
		)
end

struct SigmoidSummary
	welfare_ratio::Float64
	worst_case_welfare_ratio::Float64
	average_inference_error::Float64
	fraction_rounds_too_high_exante_allocation_prob::Float64
	average_overallocation_probability::Float64
end

function sigmoidSummaryPrint(sigmoid_summary::SigmoidSummary)
	print("Averaged welfare ratio: ")
	println(sigmoid_summary.welfare_ratio)

	print("Worst-case welfare ratio assuming that we fit p_high correctly: ")
	println(sigmoid_summary.worst_case_welfare_ratio)

	print("Average absolute inference error: ")
	println(sigmoid_summary.average_inference_error)

	print("Fraction of rounds where the ex-ante allocation probability was too high: ")
	println(sigmoid_summary.fraction_rounds_too_high_exante_allocation_prob)

	print("Expected fraction of agents who were meant to be allocated among the q fraction, but couldn't be because we made the ex-ante allocation probability too high: ")
	println(sigmoid_summary.average_overallocation_probability)
end

function plotXandY(title_name::AbstractString, x_axis_label::AbstractString, y_axis_label::AbstractString, x_axis_values::Vector{Float64}, y_axis_values::Vector{Float64})
	ourplot = plot(x_axis_values, [y_axis_values], label=[y_axis_label], lw=[1])
	plot!(ourplot, legend=:outerbottom, legendcolumns=2)
	# Note: This function is only for plotting things where the x-axis values are in [0, 1].
	xlims!(ourplot, 0, 1)
	title!(ourplot, title_name)
	xlabel!(ourplot, x_axis_label)
	ylabel!(ourplot, y_axis_label)
	println("Got to the printing phase")
	displayandpause(ourplot)
	println("^^ Actually printed")
	# TODO: This isn't actually printing a plot! Why on Earth not????
end

function plotXandYandSave(title_name::AbstractString, x_axis_label::AbstractString, y_axis_label::AbstractString, x_axis_values::Vector{Float64}, y_axis_values::Vector{Float64})
	ourplot = plot(x_axis_values, [y_axis_values], label=[y_axis_label], lw=[1])
	plot!(ourplot, legend=:outerbottom, legendcolumns=2)
	# Note: This function is only for plotting things where the x-axis values are in [0, 1].
	xlims!(ourplot, 0, 1)
	title!(ourplot, title_name)
	xlabel!(ourplot, x_axis_label)
	ylabel!(ourplot, y_axis_label)
	pathOfThisScript = @__DIR__ # Macro for acquiring directory of this script.
	file_to_be_saved = joinpath(pathOfThisScript, "..", "..", "output_images", title_name * ".svg")
	savefig(ourplot, file_to_be_saved)
end

function main()
	jokeyIntroSection()
	# testingTheSigmoidGenerator();

	# valueCDF = betaCDF(1.0, 20.0)
	valueCDF = pointMassCDF(1.0 / 21.0)
	nonzerotypes = 90
	true_cdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	# Fix parameters q and w at the start.
	# Note: We're trying to make q here a teeny bit higher than what would be at the mean for the value distribution.
	q = 1.0 / 22
	w = 0.5
	lambda = 14.0
	num_samples = 20
	num_rounds = 10
	# sigmoid_summary = sigmoidEstimationProcedure(true_cdf_values, q, w, lambda, num_samples, num_rounds)


	valueCDF = betaCDF(2.0, 10.0)
	nonzerotypes = 90
	true_cdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	# Fix parameters q and w at the start.
	q = 0.6
	w = 0.1
	lambda = 0.0
	num_samples = 20
	num_rounds = 10
	# sigmoid_summary = sigmoidEstimationProcedure(true_cdf_values, q, w, lambda, num_samples, num_rounds)
	# sigmoidSummaryPrint(sigmoid_summary)


	# Generate a graph for every lambda \in [0, +inf)
	# Iterate through w \in [0, 1]
	# // Search through all F and q values
	# // Alternative: Search through all F in Beta(alpha, beta) for bounded integer (alpha, beta)
	# Plot 1 - E_{v \sim F}[| \hat{q} - q |] that we observe from the simulation
	# // (Take the maximum E_{v \sim F}[| \hat{q} - q |] over all q and F.)

	# TODO: Currently this glitches out for lambda = +inf! Why? Shouldn't the code be robust to BR agents?
	num_samples = 10
	num_lambdas = 10
	num_ws = 20
	num_rounds = 100
	x_axis_values = Vector{Float64}(undef, num_ws)
	y_axis_values = Vector{Float64}(undef, num_ws)
	for i in 1:num_lambdas
		@assert num_lambdas > 1
		@assert num_ws > 1
		temp_x::Float64 = (i - 1) / (num_lambdas - 1)
		lambda::Float64 = Inf
		if (temp_x < 1.0)
			lambda = temp_x / (1.0 - temp_x)
		end
		for j in 1:num_ws
			w = (j - 1) / (num_ws - 1)
			# TODO: Search for the worst-case distribution F and the worst-case q.
			# valueCDF = betaCDF(2.0, 10.0)
			valueCDF = twoPointMassesCDF(0.25, 0.75, 0.5)
			nonzerotypes = 90
			true_cdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
			# Fix parameters q and w at the start.
			q = 0.6
			# @profile sigmoidEstimationProcedure(true_cdf_values, q, w, lambda, num_samples, num_rounds)
			sigmoidEstimationProcedure(true_cdf_values, q, w, lambda, num_samples, num_rounds)
			sigmoid_summary = sigmoidEstimationProcedure(true_cdf_values, q, w, lambda, num_samples, num_rounds)
			# y_axis_datum = 1.0 - sigmoid_summary.average_inference_error
			y_axis_datum = 1.0 - sigmoid_summary.average_overallocation_probability
			x_axis_values[j] = w
			y_axis_values[j] = y_axis_datum
		end
		# plotXandY("A graph", "w", "1 - E[| \\hat{q} - q |]", x_axis_values, y_axis_values)
		# plotXandY("A graph", "w", "1 - overallocation error", x_axis_values, y_axis_values)
		plotXandYandSave("Lambda = " * string(lambda) * ", q = " * string(q), "w", "1 - overallocation error", x_axis_values, y_axis_values)
	end


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
	# inferValueDistributionTest(uniformFrequency(nonzerotypes), x, lambda)
	theta_granularity = 100
	numberofsamples = 100
	x_family = straightlineDashboardFamily
	# x_family = postedpriceDashboardFamily
	# startSimulation(theta_granularity, numberofsamples, x_family, valueCDF, lambda, nonzerotypes)
	# calculateWelfareOfFamily(theta_granularity, x_family, valueCDF, lambda, nonzerotypes)

	lambda = 0.0
	nonzerotypes = 20
	theta_granularity = 10000
	x_family = straightlineDashboardFamily
	# x_family = postedpriceDashboardFamily
	# valueCDF = mixtureZeroOneCDF(0.0001, 0.0001)
	valueCDF = betaCDF(2.0, 20.0)
	# startSimulation(theta_granularity, numberofsamples, x_family, valueCDF, lambda, nonzerotypes)
	calculateWelfareOfFamily(theta_granularity, x_family, valueCDF, lambda, nonzerotypes)


	theta_granularity = 100
	starting_theta = 0.5
	q_target = 0.3
	epsilon = 0.02
	# x_family = straightlineDashboardFamily
	x_family = postedpriceDashboardFamily
	valueCDF = betaCDF(2.0, 10.0)
	true_pdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	pdfFromCDF!(true_pdf_values)
	lambda = 3.0
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	epsilon = 0.2
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	valueCDF = betaCDF(10.0, 2.0)
	true_pdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	pdfFromCDF!(true_pdf_values)
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	# theta_granularity = 100000
	# q_target = 0.00001
	# inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	# q_target = 1.0 - 0.00001
	# inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	valueCDF = betaCDF(10.0, 2.0)
	true_pdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	pdfFromCDF!(true_pdf_values)
	theta_granularity = 10000
	starting_theta = 0.001
	q_target = 0.5
	epsilon = 0.002
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	valueCDF = betaCDF(10.0, 2.0)
	true_pdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	pdfFromCDF!(true_pdf_values)
	theta_granularity = 10000
	starting_theta = 1.0 - 0.001
	q_target = 0.5
	epsilon = 0.002
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)


	valueCDF = betaCDF(2.0, 10.0)
	true_pdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	pdfFromCDF!(true_pdf_values)
	theta_granularity = 10000
	starting_theta = 0.001
	q_target = 0.5
	epsilon = 0.002
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)

	valueCDF = betaCDF(2.0, 10.0)
	true_pdf_values = valueCDF.(collect(0:nonzerotypes) ./ nonzerotypes)
	pdfFromCDF!(true_pdf_values)
	theta_granularity = 10000
	starting_theta = 1.0 - 0.001
	q_target = 0.5
	epsilon = 0.002
	inferenceQErrorWithEpsEstimationGuarantee(theta_granularity, starting_theta, q_target, epsilon, x_family, true_pdf_values, lambda)


	# x = exampleDashboard
	# valueCDF = betaCDF(2.0, 2.0)
	# lambda = 30.0
	# nonzerotypes = 90
	# displayMirage(x, valueCDF, lambda, nonzerotypes)
end

end # module MirageSimulation
