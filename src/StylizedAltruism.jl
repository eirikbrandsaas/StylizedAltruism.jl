module StylizedAltruism
using Interpolations

include("numericaltools.jl")
include("structures.jl")
include("fundamentals.jl")
include("decproblem_AHK.jl")
include("convenience.jl")

export Solve, BenchmarkParameters, switches
end # module
