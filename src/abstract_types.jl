abstract type AbstractType end

function Base.show(io::IO, data::AbstractType)
  print(io, "QuboLS.$(nameof(typeof(data)))(…)")
end

abstract type AbstractEncoding <: AbstractType end
abstract type AbstractProblem <: AbstractType end
abstract type AbstractQUBO <: AbstractType end
