using Test

#tests = ["assem", "disassem", "simulator"]
tests = ["gates", "symjit", "symbolics"]
if !isempty(ARGS)
	tests = ARGS  # Set list to same as command line args
end


dict_testfiles = Dict(
    "gates" => "Gates/gates_tests.jl",
    "symjit" => "SymJit/symjit_tests.jl",
    "symbolics" => "Symbolics/symbolics_test.jl"
)



@testset "All Tests" begin
    for t in tests
        println("Running tests: $t")
        include(dict_testfiles[t])
    end
end