module StylizedAltruism
using Interpolations

include("numericaltools.jl")
include("structures.jl")
include("fundamentals.jl")
include("decproblem_AHK.jl")
include("convenience.jl")
include("decproblem_dictator.jl")
include("decproblem_commitment.jl")

export Solve, BenchmarkParameters, switches, SolveDictator, SolveCommitment
end # module
