using PeriodicSchurBifurcationKit
using Test

@testset "PeriodicSchurBifurcationKit.jl" begin
    @testset "ODES with dense solvers" begin
		include("stuartlandau.jl")
	end

	@testset "PDEs with matrix free solvers (shooting)" begin
		include("brusselator-shooting.jl")
	end
end
