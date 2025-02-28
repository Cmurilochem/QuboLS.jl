mutable struct EncodedLinearProblem{A,B,C,D,E} <: AbstractProblem
  encoded_solution::A
  matrix::B
  rhs::C
  qubo_cost_function::D
  attributes::E
end

function encoded_linear_problem(matrix::Matrix{Float64}, rhs::Vector{Float64}, encoding::AbstractEncoding)
  (; n_variables, n_qubits, polynom) = encoding
  qubo_cost_function = Symbolics.Num(undef)

  solution = Symbolics.variables(:x, 1:n_variables)

  attributes = Dict(
    :n_qubits => n_qubits,
    :n_variables => n_variables,
    :solution => solution
  )

  return EncodedLinearProblem(
    polynom,
    matrix,
    rhs,
    qubo_cost_function,
    attributes
  )
end

function get_cost_function!(p::EncodedLinearProblem)
  (; matrix, rhs, attributes) = p
  solution = attributes[:solution]

  cost_p1 = (solution' * matrix' * matrix * solution)
  cost_p2 = -2.0 * (rhs' * matrix * solution)
  cost_p3 = rhs' * rhs

  p.attributes[:cost_function] = cost_p1 + cost_p2 + cost_p3
end

function get_qubo_cost_function!(p::EncodedLinearProblem)
  (; encoded_solution, matrix, rhs, attributes) = p
  n_variables = attributes[:n_variables]
  n_qubits = attributes[:n_qubits]

  cost_p1 = (encoded_solution' * matrix' * matrix * encoded_solution)
  cost_p1 = Symbolics.simplify(cost_p1, expand=true)

  cost_p2 = -2.0 * (rhs' * matrix * encoded_solution)
  cost_p2 = Symbolics.simplify(cost_p2, expand=true)
  s = Symbolics.variables(:s, 1:n_variables, 1:n_qubits)
  dict = Dict((s[i, j] => s[i, j]^2) for i in 1:n_variables for j in 1:n_qubits)
  cost_p2 = Symbolics.substitute(cost_p2, dict)

  cost_p3 = rhs' * rhs

  qubo_cost_function = cost_p1 + cost_p2 + cost_p3
  p.qubo_cost_function = Symbolics.simplify(qubo_cost_function, expand=true)
end

function eval_qubo_cost_function(p::EncodedLinearProblem, binary_set::Vector{Int})
  (; encoded_solution, qubo_cost_function, attributes) = p
  n_variables = attributes[:n_variables]
  n_qubits = attributes[:n_qubits]

  chopped_vectors = [collect(i) for i in IterTools.partition(binary_set, n_qubits)]
  @assert size(chopped_vectors, 1) == n_variables
  @assert size(chopped_vectors[1], 1) == n_qubits

  s = Symbolics.variables(:s, 1:n_variables, 1:n_qubits)
  dict = Dict((s[i, j] => chopped_vectors[i][j]) for i in 1:n_variables for j in 1:n_qubits)
  cost = Symbolics.substitute(qubo_cost_function, dict)
  decoded_solution = Symbolics.substitute(encoded_solution, dict)
  return Symbolics.value.(cost), Symbolics.value.(decoded_solution)
end

function eval_cost_function(p::EncodedLinearProblem, float_set::Vector{Float64})
  (; attributes) = p
  n_variables = attributes[:n_variables]
  @assert n_variables == size(float_set, 1)

  x = Symbolics.variables(:x, 1:n_variables)
  dict = Dict((x[i] => float_set[i]) for i in 1:n_variables)
  cost = Symbolics.substitute(attributes[:cost_function], dict)
  return Symbolics.value.(cost)
end
