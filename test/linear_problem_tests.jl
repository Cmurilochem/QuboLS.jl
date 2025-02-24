import QuboLS as QLS
using Test

N = 2
n_qubits = 2
range = 1.0

sol = [-3.0, 1.0]
A = [[3.0 -7.0]; [2.0 8.0]]
b = [-16.0, 2.0]

qubit_encoding = QLS.ranged_efficient_encoding(n_variables=N, n_qubits=n_qubits, range=range)
QLS.calculate_polynom!(qubit_encoding)

lp = QLS.encoded_linear_problem(A, b, qubit_encoding)
@test lp.matrix == A
@test lp.rhs == b

QLS.get_cost_function!(lp)
@test QLS.eval_cost_function(lp, sol) ≈ 0.0

QLS.get_qubo_cost_function!(lp)
cost, solution = QLS.eval_qubo_cost_function(lp, [1, 1, 1, 1])
@test cost ≈ 199.555555555
@test solution ≈ [0.6666666666666666, 0.6666666666666666]

