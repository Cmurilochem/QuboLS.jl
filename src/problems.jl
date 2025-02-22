struct LinearProblem{A,B,C} <: AbstractProblem
  solution::A
  matrix::B
  rhs::C
end

function linear_problem(solution, matrix, rhs)
  return LinearProblem(solution, matrix, rhs)
end

mutable struct EncodedLinearProblem{A,B,C,D,E,F} <: AbstractProblem
  solution::A
  encoded_solution::B
  matrix::C
  rhs::D
  qubo_cost_function::E
  cost_function::F
end

function encoded_linear_problem(matrix::Matrix{Float64}, rhs::Vector{Float64}, encoding::AbstractEncoding)
  encoded_solution = encoding.polynom
  problem_size = size(rhs, 1)
  solution = Symbolics.variables(:x, 1:problem_size)
  qubo_cost_function = Symbolics.Num(undef)
  cost_function = Symbolics.Num(undef)
  return EncodedLinearProblem(
    solution,
    encoded_solution,
    matrix,
    rhs,
    qubo_cost_function,
    cost_function
  )
end

function get_cost_function!(p::EncodedLinearProblem)
  cost_function = (p.solution' * p.matrix' * p.matrix * p.solution) - 2.0 * (p.rhs'p.matrix * p.solution) + p.rhs' * p.rhs
  p.cost_function = Symbolics.simplify(cost_function, expand=true)
end

function get_qubo_cost_function!(p::EncodedLinearProblem)
  qubo_cost_function = (p.encoded_solution' * p.matrix' * p.matrix * p.encoded_solution) - 2.0 * (p.rhs'p.matrix * p.encoded_solution) + p.rhs' * p.rhs
  p.qubo_cost_function = Symbolics.simplify(qubo_cost_function, expand=true)
end

function eval_qubo_cost_function(p::EncodedLinearProblem, binary_set::Vector{Int})
  variables = Symbolics.get_variables.(p.encoded_solution)
  n_variables = size(variables, 1)
  n_qubits = size(variables[1], 1)

  chopped_vectors = [collect(i) for i in IterTools.partition(binary_set, n_qubits)]
  @assert size(chopped_vectors, 1) == n_variables
  @assert size(chopped_vectors[1], 1) == n_qubits

  s = Symbolics.variables(:s, 1:n_variables, 1:n_qubits)
  dict = Dict((s[i, j] => chopped_vectors[i][j]) for i in 1:n_variables for j in 1:n_qubits)
  cost = Symbolics.substitute(p.qubo_cost_function, dict)
  return Symbolics.value.(cost)
end

function eval_cost_function(p::EncodedLinearProblem, float_set::Vector{Float64})
  n_variables = size(p.rhs, 1)
  @assert n_variables == size(float_set, 1)

  x = Symbolics.variables(:x, 1:n_variables)
  dict = Dict((x[i] => float_set[i]) for i in 1:n_variables)
  cost = Symbolics.substitute(p.cost_function, dict)
  return Symbolics.value.(cost)
end
