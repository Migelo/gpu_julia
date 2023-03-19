using Test
include("../src/legendre.jl")
# include("../legendre/l")

println("Testing...")

# Run the test suite
@testset "legendre_ function tests" begin
    @test legendre_(0, 0) == 1
    @test legendre_(1, 5) == 5
end
