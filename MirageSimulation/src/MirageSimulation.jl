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
	for i = reverse(2:nonzerotypes + 1)
		probabilities[i] -= probabilities[i - 1]
	end
	return probabilities
end

function mirageCDF(v_hat::Float64, x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# The value of v_hat_index is "what is the largest type index i in {0,...,n} that
	# has a b_i greater than or equal to v_hat".
	v_hat_index::Int64 = floor(v_hat * nonzerotypes)
	return mirageCDFWithIndices(v_hat_index, x, valueCDF, lambda, nonzerotypes)
end

function mirageCDFWithIndices(v_hat_index::Int64, x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# Calculates P(bid index <= `v_hat_index`).
	conditionalprobs = Vector{Float64}(undef, nonzerotypes + 1)
	for i in eachindex(conditionalprobs)
		# Compute the conditional cdf of the agent's bid, given that they have a true value index i - 1.
		conditionalprobs[i] = conditionalMirageCDFWithIndices(v_hat_index, i - 1, x, lambda, nonzerotypes)
	end
	# Calculate the probability of all true agent value regions.
	probabilities = discretizedProbabilities(valueCDF, nonzerotypes)
	# Dotting the two vectors together calculates an expectation.
	return dot(conditionalprobs, probabilities)
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

function exAnteAllocationProbability(x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)::Float64
	# This is the probability q that dashboard x allocates an item to a quantal-responding agent.
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	mirage_pdf_values = mirageCDFWithIndices.(collect(0:nonzerotypes), x, valueCDF, lambda, nonzerotypes)
	for i = reverse(2:nonzerotypes + 1)
		mirage_pdf_values[i] -= mirage_pdf_values[i - 1]
	end
	dashboard_probs = x.(agent_values)
	return dot(dashboard_probs, mirage_pdf_values)
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

function startSimulation(x::Function, valueCDF::Function, lambda::Float64, nonzerotypes::Int64)
	@assert nonzerotypes >= 1
	@assert lambda >= 0.0
	println("Simulation Started!")
	agent_values = collect(0:nonzerotypes) ./ nonzerotypes
	true_cdf_values = valueCDF.(agent_values)
	mirage_cdf_values = mirageCDFWithIndices.(collect(0:nonzerotypes), x, valueCDF, lambda, nonzerotypes)
	print("Ex ante allocation probability (q): ")
	println(exAnteAllocationProbability(x, valueCDF, lambda, nonzerotypes))
	# If you were to instead use the command
	# mirage_cdf_values = mirageCDF.(agent_values, x, valueCDF, lambda, nonzerotypes)
	# it sometimes results in a floating point rounding error that causes some adjacent v_hat values to
	# produce the same output. And that's bad!
	plotResults(agent_values, true_cdf_values, mirage_cdf_values)
end

function main()
	jokeyIntroSection()
	x = identityDashboard
	valueCDF = betaCDF(2.0, 2.0)
	lambda = 30.0
	nonzerotypes = 90
	startSimulation(x, valueCDF, lambda, nonzerotypes)
	x = exampleDashboard
	valueCDF = betaCDF(2.0, 2.0)
	lambda = 30.0
	nonzerotypes = 90
	startSimulation(x, valueCDF, lambda, nonzerotypes)
end

end # module MirageSimulation
