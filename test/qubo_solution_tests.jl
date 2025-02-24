using QuboLS
using Random
using JuMP, QUBO
using Test

Random.seed!(42)

n_variables = 3
n_qubits = 12
range = 3.0

A = rand(n_variables, n_variables)
b = rand(n_variables)
exact_solution = A \ b

@test isapprox(A * exact_solution, b)

encoding = QuboLS.ranged_efficient_encoding(
  n_variables=n_variables, n_qubits=n_qubits, range=range
)
QuboLS.calculate_polynom!(encoding)

lp = QuboLS.encoded_linear_problem(A, b, encoding)
QuboLS.get_qubo_cost_function!(lp)

qubo = QuboLS.QUBO(lp)
QuboLS.get_qubo_matrix!(qubo)

model = JuMP.Model(() -> ToQUBO.Optimizer(ExactSampler.Optimizer))
JuMP.@variable(model, s[1:n_variables*n_qubits], Bin)
JuMP.@objective(model, Min, s' * qubo.matrix * s + qubo.offset)

JuMP.optimize!(model)

qubo_binary_solution = Int.(JuMP.value.(s))
qubo_energy = JuMP.objective_value(model)
my_cost, approx_solution = QuboLS.eval_qubo_cost_function(lp, qubo_binary_solution)

@test isapprox(qubo_energy, my_cost)
@test isapprox(approx_solution, exact_solution; atol=0.005)
