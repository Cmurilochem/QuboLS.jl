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

function get_qubo_matrix!(q::QUBO; check_matrix::Bool=false)
  n_variables = q.attributes[:n_variables]
  n_qubits = q.attributes[:n_qubits]
  cost = q.attributes[:qubo_cost_function]

  s = Symbolics.variables(:s, 1:n_variables, 1:n_qubits)
  q_matrix = q.matrix

  for n in 1:(n_variables * n_qubits)
      i, j = divrem(n - 1, n_qubits) .+ 1
      for m in n:(n_variables * n_qubits)
          k, l = divrem(m - 1, n_qubits) .+ 1
          term = s[i, j] * s[k, l]
          coeff = Symbolics.coeff(cost, term)
          if !isapprox(coeff, 0)
              q_matrix[n, m] = coeff
          end
      end
  end

  #for (i, j) in Iterators.product(1:n_variables, 1:n_qubits)
  #  for (k, l) in Iterators.product(1:n_variables, 1:n_qubits)
  #    n = (i - 1) * n_qubits + j
  #    m = (k - 1) * n_qubits + l
  #    term = s[i, j] * s[k, l]
  #    coeff = Symbolics.coeff(cost, term)
  #    if m >= n && !isapprox(coeff, 0)
  #      q_matrix[n, m] = coeff
  #    end
  #  end
  #end

  q.offset = Symbolics.coeff(cost)

  if check_matrix
    s_total = Vector{Symbolics.Num}(undef, 0)
    for i in 1:n_variables
      term = [s[i, j] for j in 1:n_qubits]
      s_total = vcat(s_total, term)
    end
    calc_cost_from_q = s_total' * q.matrix * s_total + q.offset
    calc_cost_from_q = Symbolics.simplify(calc_cost_from_q, expand=true)

    @assert calc_cost_from_q - cost == 0
  end

end
