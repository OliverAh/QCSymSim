import Symbolics

@testset "Complex Symbols" begin
    println("Running complex symbols tests...")
    
    @testset "trig of complex symbol" begin
        Symbolics.@variables y x::Complex
        # typeof(im*y)
        # typeof(x)
        sin(im*y)
        @test_throws TypeError sin(x)
    end
end
