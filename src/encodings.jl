mutable struct RangedEfficientEncoding{A,B,C,D,E} <: AbstractEncoding
  n_variables::A
  n_qubits::B
  range::C
  var_base_name::D
  polynom::E
end

function ranged_efficient_encoding(;
  n_variables::Int=1,
  n_qubits::Int=2,
  range::Union{Vector{Float64},Float64}=1.0,
  var_base_name::Symbol=:s
)
  if !isa(range, Vector)
    range = fill(range, n_qubits)
  end
  polynom = Vector{Symbolics.Num}(undef, n_variables)
  return RangedEfficientEncoding(n_variables, n_qubits, range, var_base_name, polynom)
end

function calculate_polynom!(encoding::RangedEfficientEncoding)
  (; n_variables, n_qubits, range, var_base_name) = encoding

  int_max = 2^n_qubits - 1

  @assert isa(range, Vector)
  max_absval = [r / int_max for r in range]

  var_set = Symbolics.variables(var_base_name, 1:n_variables, 1:n_qubits)
  coeffs = Vector{Float64}(undef, n_qubits)

  coeffs[1] = -2^(n_qubits - 1) * max_absval[1]
  coeffs[2:end] = 2 .^ (2:n_qubits) .* max_absval[2:end]

  encoding.polynom = var_set * coeffs
end
