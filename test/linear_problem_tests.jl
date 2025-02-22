using QuboEncoding
using Symbolics
using Test

N = 2
n_qubits = 2
range = 1.0

sol = [-3.0, 1.0]
A = [[3.0 -7.0]; [2.0 8.0]]
b = [-16.0, 2.0]

lp = QuboEncoding.linear_problem(sol, A, b)
@test lp isa QuboEncoding.LinearProblem

qubit_encoding = QuboEncoding.ranged_efficient_encoding(n_variables=N, n_qubits=n_qubits, range=range)
QuboEncoding.calculate_polynom!(qubit_encoding)

lp = QuboEncoding.encoded_linear_problem(A, b, qubit_encoding)

QuboEncoding.get_cost_function!(lp)
@test QuboEncoding.eval_cost_function(lp, sol) ≈ 0.0 

QuboEncoding.get_qubo_cost_function!(lp)
@test QuboEncoding.eval_qubo_cost_function(lp, [1, 1, 1, 1]) ≈ 199.555555555
