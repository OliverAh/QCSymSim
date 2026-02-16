import Symbolics
import ..BitsRegs.AbstractBit
#import ..BitsRegs.Bit
import ..BitsRegs.MapBitID
import ..BitsRegs.QBit

abstract type Bit <: AbstractBit end

struct _00_Gate <: AbstractInternalSingleQubitQuantumGate
    function _00_Gate(;name::AbstractString, qubits_t::Array{Int,1}, step::UInt, is_treat_numeric_only::Bool)
        matrix_numeric = [1 0; 0 0]
        base_gate = make_BaseQuantumGate(name=name, name_short="_00", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
       return base_gate
    end
end

struct _11_Gate <: AbstractInternalSingleQubitQuantumGate
    function _11_Gate(;name::AbstractString, qubits_t::Array{Int,1}, step::UInt, is_treat_numeric_only::Bool)
        matrix_numeric = [0 0; 0 1]
        base_gate = make_BaseQuantumGate(name=name, name_short="_11", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
       return base_gate
    end
end

struct I_Gate <: AbstractSingleQubitQuantumGate
    function I_Gate(;name::AbstractString, qubits_t::Array{Int,1}, step::UInt, is_treat_numeric_only::Bool)
        matrix_numeric = [1 0; 0 1]
        base_gate = make_BaseQuantumGate(name=name, name_short="I", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
       return base_gate
    end
end

struct X_Gate <: AbstractSingleQubitQuantumGate
    function X_Gate(;name::AbstractString, qubits_t::Array{Int,1}, step::UInt, is_treat_numeric_only::Bool)
        matrix_numeric = [0 1; 1 0]
        base_gate = make_BaseQuantumGate(name=name, name_short="X", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
        return base_gate
    end
end

struct H_Gate{T<:AbstractBit} <: AbstractSingleQubitQuantumGate
    function H_Gate{T}(;name::AbstractString, qubits_t::AbstractArray{T,1}, step::UInt, is_treat_numeric_only::Bool) where {T<:AbstractBit}
        matrix_numeric = (1/sqrt(2))*[1 1; 1 -1]
        base_gate = make_BaseQuantumGate(name=name, name_short="H", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
        base_gate.matrix_alt = 1/Symbolics.ssqrt(2)*Matrix{Symbolics.Num}([1 1; 1 -1])
        #base_gate.matrix_alt = Symbolics.scalarize(base_gate.matrix_alt)
        return base_gate
    end
end
function H_Gate_outer(;name::AbstractString, qubits_t::Array{Int,1}, step::UInt, is_treat_numeric_only::Bool)
    println("ended up in OUTER constructor")
    context = MapBitID()
    qubits_t_new = collect(QBit(context=context) for q in qubits_t)
    println(typeof(qubits_t_new))
    return H_Gate{typeof(qubits_t_new[1])}(name=name, qubits_t=qubits_t_new, step=step, is_treat_numeric_only=is_treat_numeric_only)
end

struct CX_Gate <: AbstractMultiQubitQuantumGate
    function CX_Gate(;name::AbstractString, qubits_t::Array{Int,1}, qubits_c::Array{Int,1}, step::UInt, is_treat_numeric_only::Bool)
        gates_t = Dict(qubits_t[1] => (I_Gate, X_Gate))
        gates_c = Dict(qubits_c[1] => (_00_Gate, _11_Gate))
        println(nameof(I_Gate))
        
        base_gate = make_MultiQubitQuantumGate(name=name, name_short="CX",qubits_t=qubits_t, qubits_c=qubits_c, step=step, is_treat_numeric_only=is_treat_numeric_only,
        parameters=nothing, gates_t=gates_t, gates_c=gates_c)
        return base_gate
    end
end