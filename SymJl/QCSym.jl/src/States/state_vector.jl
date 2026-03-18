import Symbolics
abstract type AbstractState end
abstract type AbstractQuantumState <: AbstractState end

abstract type AbstractStateVector <: AbstractQuantumState end

struct StateVector{T} <: AbstractStateVector
    vector_symbolic::AbstractVector{<:Symbolics.Num}
    vector_numeric::AbstractVector{T}
end

function StateVector(T::Type, num_qubits::Int)
        T = T === Nothing ? ComplexF64 : T
        _vector_numeric = zeros(T, 2^num_qubits)
        _vector_numeric[1] = 1.0
        _vector_symbolic = Symbolics.scalarize(Symbolics.@variables(φ[1:2^num_qubits])[1])
        StateVector{T}(_vector_symbolic, _vector_numeric)
end


function StateVector(T::Type, kets::AbstractVector{TT}) where {TT<:Bool}
    T = T === Nothing ? ComplexF64 : T
    ket_0 = Array{T}([1.0, 0.0])
    ket_1 = Array{T}([0.0, 1.0])
    _vector_numeric = kets[1] ? ket_1 : ket_0
    for i in 2:length(kets)
        _vector_numeric = kron(_vector_numeric, kets[i] ? ket_1 : ket_0)
    end
    _vector_symbolic = Symbolics.scalarize(Symbolics.@variables(φ[1:2^num_qubits])[1])
    StateVector{T}(_vector_symbolic, _vector_numeric)
end

function StateVector(vec::Vector{T}) where {T<:Number}
        _vector_numeric = vec
        num_qubits = round(Int, log2(length(vec)))
        _vector_symbolic = Symbolics.scalarize(Symbolics.@variables(φ[1:2^num_qubits])[1])
        StateVector{T}(_vector_symbolic, _vector_numeric)
end

function StateVector(kets::Vector{Vector{TT}}) where {TT<:Complex}
        _vector_numeric = kets[1]
        num_qubits = round(Int, log2(length(kets)))
        for i in 2:length(kets)
            _vector_numeric = kron(_vector_numeric, kets[i])
        end
        _vector_symbolic = Symbolics.scalarize(Symbolics.@variables(φ[1:2^num_qubits])[1])
        StateVector{TT}(_vector_symbolic, _vector_numeric)
end