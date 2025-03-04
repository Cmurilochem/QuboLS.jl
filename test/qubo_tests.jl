import QuboLS as QLS
import LinearAlgebra
using Test

N = 2
n_qubits = 4
range = 1.0

sol = [-3.0, 1.0]
A = [[3.0 -7.0]; [2.0 8.0]]
b = [-16.0, 2.0]

qubit_encoding = QLS.ranged_efficient_encoding(n_variables=N, n_qubits=n_qubits, range=range)
QLS.calculate_polynom!(qubit_encoding)

lp = QLS.encoded_linear_problem(A, b, qubit_encoding)
QLS.get_qubo_cost_function!(lp)

qubo = QLS.QUBO(lp)
QLS.get_qubo_matrix!(qubo, check_matrix=true)

@test qubo.offset == 260.0
@test LinearAlgebra.istriu(qubo.matrix)
