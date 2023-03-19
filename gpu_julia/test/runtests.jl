using Test
include("../src/legendre.jl")

println("Testing...")

@testset "legendre_ function tests" begin
    @test legendre_(0, 0) == 1
    @test legendre_(1, 5) == 5
    @test legendre_(2, 0.5) == 0.5 * (3 * 0.5 * 0.5 - 1)
    @test legendre_(3, 0.2) == 0.5 * (5 * 0.2^3 - 3 * 0.2)
    @test legendre_(4, -0.2) == (1.0 / 8) * (35 * (-0.2)^4 - 30 * (-0.2)^2 + 3)
    @test legendre_(5, -0.89) == (1.0 / 8) * (63 * (-0.89)^5 - 70 * (-0.89)^3 + 15 * (-0.89))
    @test legendre_(6, 0.24) == (1.0 / 16) * (
        231 * 0.24^6 - 315 * 0.24^4 + 105 * 0.24^2 - 5
    )
    @test legendre_(7, 0.54) == (1.0 / 16) * (
        429 * 0.54^7 - 693 * 0.54^5 + 315 * 0.54^3 - 35 * 0.54
    )
    @test legendre_(8, 0.66) == (1.0 / 128) * (6435 * 0.66^8
                                               -
                                               12012 * 0.66^6
                                               +
                                               6930 * 0.66^4
                                               -
                                               1260 * 0.66^2
                                               +
                                               35)
    @test legendre_(9, 0.31) == (1.0 / 128) * (
        12155 * 0.31^9
        -
        25740 * 0.31^7
        +
        18018 * 0.31^5
        -
        4620 * 0.31^3
        +
        315 * 0.31
    )
    @test legendre_(10, 0.12) == (1.0 / 256) * (
        46189 * 0.12^10
        -
        109395 * 0.12^8
        +
        90090 * 0.12^6
        -
        30030 * 0.12^4
        +
        3465 * 0.12^2
        -
        63
    )
    @test_throws AssertionError legendre_(-1, 5)
    @test_throws AssertionError legendre_(11, 5)
end

@testset "legendre" begin
    @test sqrt(2 * 0 + 1) * legendre_(0, 0.5) == legendre(0, 0.5)
    @test sqrt(2 * 1 + 1) * legendre_(1, 0.5) == legendre(1, 0.5)
    @test sqrt(2 * 2 + 1) * legendre_(2, 0.5) == legendre(2, 0.5)
    @test sqrt(2 * 3 + 1) * legendre_(3, 0.5) == legendre(3, 0.5)
    @test sqrt(2 * 4 + 1) * legendre_(4, 0.5) == legendre(4, 0.5)
    @test sqrt(2 * 5 + 1) * legendre_(5, 0.5) == legendre(5, 0.5)
    @test sqrt(2 * 6 + 1) * legendre_(6, 0.5) == legendre(6, 0.5)
    @test sqrt(2 * 7 + 1) * legendre_(7, 0.5) == legendre(7, 0.5)
    @test sqrt(2 * 8 + 1) * legendre_(8, 0.5) == legendre(8, 0.5)
    @test sqrt(2 * 9 + 1) * legendre_(9, 0.5) == legendre(9, 0.5)
    @test sqrt(2 * 10 + 1) * legendre_(10, 0.5) == legendre(10, 0.5)

end
