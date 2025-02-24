import QuboLS as QLS
import Random
import JuMP, QUBO, ToQUBO, QUBODrivers
using Test

Random.seed!(42)

# Define different parameter sets to test
test_cases = [
  (n_variables=2, n_qubits=6, range=2.0),
  (n_variables=3, n_qubits=6, range=10.0),
  (n_variables=4, n_qubits=6, range=7.0),
]

for params in test_cases
  @testset "Testing with n_variables=$(params.n_variables), n_qubits=$(params.n_qubits), range=$(params.range)" begin

    # Extract parameters
    n_variables = params.n_variables
    n_qubits = params.n_qubits
    range = params.range

    # Generate linear system
    A = rand(n_variables, n_variables)
    b = rand(n_variables)
    exact_solution = A \ b

    # Test if exact solution satisfies Ax = b
    @test isapprox(A * exact_solution, b)

    # Encode problem using QUBO
    encoding = QLS.ranged_efficient_encoding(
      n_variables=n_variables, n_qubits=n_qubits, range=range
    )
    QLS.calculate_polynom!(encoding)

    lp = QLS.encoded_linear_problem(A, b, encoding)
    QLS.get_qubo_cost_function!(lp)

    qubo = QLS.QUBO(lp)
    QLS.get_qubo_matrix!(qubo)

    # Solve QUBO using JuMP and QUBO.jl
    sampler = QUBODrivers.ExactSampler.Optimizer
    optimizer = ToQUBO.Optimizer(sampler)

    model = JuMP.Model(() -> optimizer)
    JuMP.@variable(model, s[1:n_variables*n_qubits], Bin)
    JuMP.@objective(model, Min, s' * qubo.matrix * s + qubo.offset)

    JuMP.optimize!(model)

    # Retrieve solution
    qubo_binary_solution = Int.(JuMP.value.(s))
    qubo_energy = JuMP.objective_value(model)
    my_cost, approx_solution = QLS.eval_qubo_cost_function(lp, qubo_binary_solution)

    # Validate QUBO solution
    @test isapprox(qubo_energy, my_cost)
    @test isapprox(approx_solution, exact_solution; atol=1)
  end
end
