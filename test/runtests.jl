module QuboEncodingTest

using Test

@testset "QuboEncoding" begin
  @testset "encoding" begin
    include("encodings_tests.jl")
  end
  @testset "linear problem" begin
    include("linear_problem_tests.jl")
  end
end

end
