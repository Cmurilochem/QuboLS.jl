mutable struct QUBO{A,B,C} <: AbstractQUBO
  matrix::A
  offset::B
  attributes::C
end

function QUBO(p::EncodedLinearProblem)
  n_variables = p.attributes[:n_variables]
  n_qubits = p.attributes[:n_qubits]
  matrix = zeros(n_qubits * n_variables, n_qubits * n_variables)

  qubo_cost_function = p.qubo_cost_function
  offset = 0.0
  encoded_solution = p.encoded_solution

  attributes = Dict(
    :n_qubits => n_qubits,
    :n_variables => n_variables,
    :qubo_cost_function => qubo_cost_function,
    :encoded_solution => encoded_solution,
  )

  return QUBO(matrix, offset, attributes)
end

function get_qubo_matrix!(q::QUBO)
  n_variables = q.attributes[:n_variables]
  n_qubits = q.attributes[:n_qubits]
  cost = q.attributes[:qubo_cost_function]

  s = Symbolics.variables(:s, 1:n_variables, 1:n_qubits)

  for (i, j) in Iterators.product(1:n_variables, 1:n_qubits)
    for (k, l) in Iterators.product(1:n_variables, 1:n_qubits)
      n = (i - 1) * n_qubits + j
      m = (k - 1) * n_qubits + l
      term = s[i, j] * s[k, l]
      coeff = Symbolics.coeff(cost, term)
      if m >= n && !isapprox(coeff, 0)
        q.matrix[n, m] = coeff
      end
    end
  end

  q.offset = Symbolics.coeff(cost)
end
