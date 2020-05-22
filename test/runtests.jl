using StylizedAltruism
using Test

include("test_nesting.jl")

@testset "StylizedAltruism.jl" begin
    # Write your own tests here.
    # Redundant test just to check whether the package development worked (deleted later)
   @test StylizedAltruism.ModPar().β > 0
end

@testset "Does changing switches work as expected?" begin
    @test test_altr_vs_eta() < 0.005 # Allow a small deviation due to numerical issues
end
