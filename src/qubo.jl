mutable struct QUBO{A,B,C} <: AbstractQUBO
  matrix::A
  offset::B
  attributes::C
end

function QUBO(p::EncodedLinearProblem)
  (; qubo_cost_function, encoded_solution, attributes) = p
  n_variables = attributes[:n_variables]
  n_qubits = attributes[:n_qubits]

  matrix = zeros(n_qubits * n_variables, n_qubits * n_variables)
  offset = 0.0

  attributes = Dict(
    :n_qubits => n_qubits,
    :n_variables => n_variables,
    :qubo_cost_function => qubo_cost_function,
    :encoded_solution => encoded_solution,
  )
  return QUBO(matrix, offset, attributes)
end

function get_qubo_matrix!(q::QUBO; check_matrix::Bool=false)
  (; matrix, offset, attributes) = q
  n_variables = attributes[:n_variables]
  n_qubits = attributes[:n_qubits]
  cost = attributes[:qubo_cost_function]

  s = Symbolics.variables(:s, 1:n_variables, 1:n_qubits)

  @inbounds for (i, j) in Iterators.product(1:n_variables, 1:n_qubits)
    n = (i - 1) * n_qubits + j
    s_ij = s[i, j]

    for (k, l) in Iterators.product(1:n_variables, 1:n_qubits)
      m = (k - 1) * n_qubits + l
      term = s_ij * s[k, l]

      if m >= n
        coeff = Symbolics.coeff(cost, term)
        @assert !isapprox(coeff, 0)
        matrix[n, m] = coeff
      end

    end
  end

  q.matrix = matrix
  offset = Symbolics.coeff(cost)
  q.offset = offset

  if check_matrix
    s_total = Vector{Symbolics.Num}(undef, 0)
    for i in 1:n_variables
      term = [s[i, j] for j in 1:n_qubits]
      s_total = vcat(s_total, term)
    end

    calc_cost_from_q = s_total' * matrix * s_total + offset
    calc_cost_from_q = Symbolics.simplify(calc_cost_from_q, expand=true)

    @assert calc_cost_from_q - cost == 0
  end

end
