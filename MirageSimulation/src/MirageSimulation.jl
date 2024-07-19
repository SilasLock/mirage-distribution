module MirageSimulation
# This is the main simulation script. It's meant to be called from "init.jl" and nowhere else.
# Note that in order for Julia's package manager to run the project properly, this file needs
# to be named "MergeSimulation.jl". Otherwise, Julia claims the package/project isn't installed. 

# Grab solely the "Beta" and "cdf" functions from the Distributions package, and nothing else.
# Doing so avoids polluting the namespace.
using Distributions: Beta, cdf

function kumaraswamyCDF(x::Float64, a::Float64, b::Float64)
	@assert a > 0
	@assert b > 0
	@assert x >= 0
	@assert x <= 1
	return 1.0 - ((1.0 - (x^a))^(b - 1.0))
end

function betaCDF(x::Float64, alpha::Float64, beta::Float64)
    return cdf(Beta(alpha, beta), x)
end

function startSimulation()
	println("Simulation code will need to be placed here.")
	println("For now, here's two outputs of our CDFs with a totally arbitrary parameterization I gave them.")
	print("Kumaraswamy distribution: ")
	println(kumaraswamyCDF(0.5, 3.0, 6.0))
	print("Beta distribution: ")
	println(betaCDF(0.6, 2.0, 3.0))
	print("Are these numbers probabilities between 0 and 1? Then great!")
	println(" We've already accomplished the first step of the project.")
end


function main()
	println("Simulator Started!")
	startSimulation()
end

end # module MirageSimulation
