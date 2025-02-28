import QuboLS as QLS
import JuMP, QUBO
using Test

A = [[3.0 -7.0]; [2.0 8.0]]
b = [-16.0, 2.0]
exact_solution = [-3.0, 1.0]

n_variables = size(b, 1)
n_qubits = 8
range = 6.5

encoding = QLS.ranged_efficient_encoding(
  n_variables=n_variables, n_qubits=n_qubits, range=range
)

optimizer = QUBO.ExactSampler.Optimizer
model = JuMP.Model(optimizer)

solution, objective_function = QLS.solve(encoding, A, b, model)

@test isapprox(solution, exact_solution; atol=0.1)
@test isapprox(objective_function, 0.0; atol=0.1)

