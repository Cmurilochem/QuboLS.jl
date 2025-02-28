module QuboLSTest
using Test

@testset "binary encoding" begin
  include("encodings_tests.jl")
end
@testset "linear problem encoding" begin
  include("linear_problem_tests.jl")
end
@testset "qubo matrix set up" begin
  include("qubo_tests.jl")
end
@testset "qubo solution" begin
  include("qubo_solution_tests.jl")
end
@testset "qubo interface" begin
  include("interfaces_test.jl")
end
end
