function solve(
  encoding::AbstractEncoding,
  A::Matrix{Float64},
  b::Vector{Float64},
  model::JuMP.Model
)
  (; n_variables, n_qubits) = encoding
  QLS.calculate_polynom!(encoding)

  encoded_linear_problem = QLS.encoded_linear_problem(A, b, encoding)
  QLS.get_qubo_cost_function!(encoded_linear_problem)

  qubo_problem = QLS.QUBO(encoded_linear_problem)
  QLS.get_qubo_matrix!(qubo_problem)

  JuMP.@variable(model, s[1:n_variables*n_qubits], Bin)
  JuMP.@objective(model, Min, s' * qubo_problem.matrix * s + qubo_problem.offset)
  JuMP.optimize!(model)

  qubo_binary_solution = Int.(JuMP.value.(s))
  obj_value = JuMP.objective_value(model)
  cost, approx_solution = QLS.eval_qubo_cost_function(encoded_linear_problem, qubo_binary_solution)

  @assert isapprox(obj_value, cost)
  return approx_solution, obj_value
end

