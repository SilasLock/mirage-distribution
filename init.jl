# This is the initialization script. Its sole purpose is to call main() from the MirageSimulation.jl file.
# To do so, navigate to the directory containing this init.jl file and then use the command
# julia --project=MirageSimulation init.jl
# to run the program.

# This activates the MirageSimulation environment, which must be located in the same
# directory as this init.jl script.
using Pkg: activate, instantiate
pathOfThisScript = @__DIR__ # Macro for acquiring directory of this script.
activate(pathOfThisScript * "\\" * "MirageSimulation", io=devnull)
# Observe that we suppress the annoying text usually caused with activate() by using io=devnull.
instantiate()
# We're, at present, not choosing to suppress the precompilation progress bar that
# comes with the instantiate() function.

using MirageSimulation
MirageSimulation.main()
