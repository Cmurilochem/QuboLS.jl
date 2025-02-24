module QuboLS
const QLS = QuboLS

using LinearAlgebra
using Symbolics
using IterTools

export ranged_efficient_encoding
export calculate_polynom!
export encoded_linear_problem
export get_cost_function!
export eval_cost_function
export get_qubo_cost_function!
export eval_qubo_cost_function
export QUBO
export get_qubo_matrix!

include("abstract_types.jl")
include("encodings.jl")
include("problems.jl")
include("qubo.jl")

end
