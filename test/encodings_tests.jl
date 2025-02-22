using QuboEncoding
using LinearAlgebra
using Symbolics
using Test

n_variables = 1
n_qubits = 2
range = 1.0

encoding = QuboEncoding.ranged_efficient_encoding(n_variables=n_variables, n_qubits=n_qubits, range=range)
@test encoding isa QuboEncoding.RangedEfficientEncoding
@test encoding.n_qubits == n_qubits
@test encoding.range == fill(range, n_qubits)
@test encoding.var_base_name == :s

QuboEncoding.calculate_polynom!(encoding)
@test encoding.polynom isa Vector{Symbolics.Num}

s = Symbolics.get_variables(encoding.polynom[1])
V = Symbolics.substitute(encoding.polynom[1], Dict(s[1] => 1, s[2] => 1))
@test Symbolics.value.(V) ≈ 0.666666666

n_variables = 2

s = Symbolics.variables(:s, 1:n_variables, 1:n_qubits)

encoding = QuboEncoding.ranged_efficient_encoding(n_variables=n_variables)
@test encoding isa QuboEncoding.RangedEfficientEncoding
@test encoding.n_qubits == n_qubits
@test encoding.range == fill(range, n_qubits)
@test encoding.var_base_name == :s

QuboEncoding.calculate_polynom!(encoding)
@test encoding.polynom isa Vector{Symbolics.Num}

V = Symbolics.substitute(encoding.polynom, Dict(s[1, 1] => 1, s[1, 2] => 1, s[2, 1] => 1, s[2, 2] => 1))
@test Symbolics.value.(V) ≈ [0.666666666, 0.666666666]
