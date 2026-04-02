using Test

import Symbolics
import SymbolicUtils
import SymJit

@testset "Complex Symbols" begin
    println("Running Complex Symbols tests...")
    
    @testset "x::Complex" begin
        @testset "scalar real" begin
            Symbolics.@variables x::Complex

            expr = x^2 + 2x + 1
            f = SymJit.compile_func([x], [expr]; dtype=:complex)

            x_num = 2.0
            expected = 9.0+0.0im
            as_is = f([x_num])
            @test length(as_is) == 1 && real(as_is[1]) ≈ real(expected) && imag(as_is[1]) == imag(expected)
        end

        @testset "scalar imag" begin
            Symbolics.@variables x::Complex

            expr = x^2 + 2x + 1
            f = SymJit.compile_func([x], [expr]; dtype=:complex)

            x_num = 2.0 + 3.0im
            expected = 0.0 + 18.0im
            as_is = f([x_num])
            @test length(as_is) == 1 && real(as_is[1]) == 0.0 && imag(as_is[1]) == 18.0
        end

        @testset "scalar trig real" begin
            Symbolics.@variables x

            expr = sin(x)
            f = SymJit.compile_func([x], [expr]; dtype=:complex)

            x_num = 2.0
            expected = sin(2.0)*cosh(0.0) + im*cos(2.0)*sinh(0.0)
            as_is = f([x_num])
            @test real(as_is[1]) ≈ real(expected) && imag(as_is[1]) == imag(expected) && length(as_is) == 1
        end

        @testset "scalar trig imag" begin
            Symbolics.@variables x::Complex

            @test_throws TypeError sin(x)
        end

        @testset "vector scalarized real" begin
            x = Symbolics.@variables(x[1:2])[1]

            xx = Symbolics.scalarize(x)
            expr = xx[1]*xx[2] + xx[1]
            f = SymJit.compile_func([xx[1], xx[2]], [expr]; dtype=:complex)
            xx_num = [2.0 , 3.0]
            expected = 8.0+0.0im
            as_is = f([xx_num[1], xx_num[2]])
            @test length(as_is) == 1 && real(as_is[1]) ≈ real(expected) && imag(as_is[1]) == imag(expected)
        end

        @testset "vector scalarized imag" begin
            x = Symbolics.@variables((x::Complex)[1:2])[1]

            xx = Symbolics.scalarize(x)
            expr = xx[1]*xx[2] + xx[1]
            f = SymJit.compile_func([xx[1], xx[2]], [expr]; dtype=:complex)
            xx_num = [2.0+1.0im, 3.0+2im]
            expected = 6.0+8.0im
            as_is = f([xx_num[1], xx_num[2]])
            @test length(as_is) == 1 && real(as_is[1]) ≈ real(expected) && imag(as_is[1]) == imag(expected)
        end

        
    end
end

