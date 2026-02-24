import Symbolics
import ..BitsRegs.AbstractBit
#import ..BitsRegs.Bit
import ..BitsRegs.MapBitID
import ..BitsRegs.QBit

struct _00_Gate{T<:AbstractBit} <: AbstractInternalSingleQubitQuantumGate{T}
    
    function _00_Gate(;name::AbstractString, qubits_t::AbstractVector{<:AbstractBit}, step::Int, is_treat_numeric_only::Bool)
        matrix_numeric = [1 0; 0 0]
        base_gate = make_BaseQuantumGate(name=name, name_short="_00", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
       return base_gate
    end
end

struct _11_Gate{T<:AbstractBit} <: AbstractInternalSingleQubitQuantumGate{T}
    function _11_Gate(;name::AbstractString, qubits_t::AbstractVector{<:AbstractBit}, step::Int, is_treat_numeric_only::Bool)
        matrix_numeric = [0 0; 0 1]
        base_gate = make_BaseQuantumGate(name=name, name_short="_11", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
       return base_gate
    end
end

struct I_Gate{T<:AbstractBit} <: AbstractSingleQubitQuantumGate{T}
    function I_Gate(;name::AbstractString, qubits_t::AbstractVector{<:AbstractBit}, step::Int, is_treat_numeric_only::Bool)
        matrix_numeric = [1 0; 0 1]
        base_gate = make_BaseQuantumGate(name=name, name_short="I", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
       return base_gate
    end
end

struct X_Gate{T<:AbstractBit} <: AbstractSingleQubitQuantumGate{T}
    function X_Gate(;name::AbstractString, qubits_t::AbstractVector{<:AbstractBit}, step::Int, is_treat_numeric_only::Bool)
        matrix_numeric = [0 1; 1 0]
        base_gate = make_BaseQuantumGate(name=name, name_short="X", shape=(2,2),
        qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
        parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
        ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
        return base_gate
    end
end

struct H_Gate{T<:AbstractBit} <: AbstractSingleQubitQuantumGate{T}
    @insert_fields_AbstractQuantumGate()
    @constructor_from_mutable_base(H_Gate, mutable_BaseQuantumGate_for_construction)
    
    # function H_Gate{T}(;name::AbstractString, qubits_t::AbstractVector{T}, step::Int, is_treat_numeric_only::Bool) where {T<:AbstractBit}
    #     matrix_numeric = (1/sqrt(2))*[1 1; 1 -1]
    #     base_gate = make_BaseQuantumGate(name=name, name_short="H", shape=(2,2),
    #     qubits_t=qubits_t, qubits_c=nothing, step=step, num_summands_decomposed=1,
    #     parameters=nothing, is_treat_numeric_only=is_treat_numeric_only,
    #     ids_matrix_zeros=nothing, matrix_numeric=matrix_numeric)
    #     base_gate.matrix_alt = 1/Symbolics.ssqrt(2)*Matrix{Symbolics.Num}([1 1; 1 -1])
    #     #base_gate.matrix_alt = Symbolics.scalarize(base_gate.matrix_alt)
    #     return base_gate
    # end
end

function H_Gate_for_Circuit(;name_prefix::String="", qubits_t::AbstractVector{QBit}, step::Int, is_treat_numeric_only::Bool)
    println("construct H_Gate for Circuit")
    base_gate = mutable_BaseQuantumGate_for_construction(is_parametric=false, is_treat_numeric_only=is_treat_numeric_only, 
        name_prefix=name_prefix, name_short="H", qubits_t=qubits_t, qubits_c=nothing,
        step=step, num_summands_decomposed=1)
    
     return H_Gate(base_gate)
end

struct CX_Gate{T<:AbstractBit} <: AbstractMultiQubitQuantumGate{T}
    @insert_fields_AbstractQuantumGate()
    @constructor_from_mutable_base(CX_Gate, mutable_BaseQuantumGate_for_construction)
end

function CX_Gate_for_Circuit(;name_prefix::String="", qubits_t::AbstractVector{QBit}, qubits_c::AbstractVector{QBit}, step::Int, is_treat_numeric_only::Bool)
    println("construct CX_Gate for Circuit")
    base_gate = mutable_BaseQuantumGate_for_construction(is_parametric=false, is_treat_numeric_only=is_treat_numeric_only, 
        name_prefix=name_prefix, name_short="CX", qubits_t=qubits_t, qubits_c=qubits_c,
        step=step, num_summands_decomposed=2)
    base_gate.matrix_numeric = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
    
     return CX_Gate(base_gate)
end


# struct CX_Gate{T<:AbstractBit} <: AbstractMultiQubitQuantumGate{T}
#     function CX_Gate(;name::AbstractString, qubits_t::AbstractVector{AbstractBit}, qubits_c::AbstractVector{AbstractBit}, step::Int, is_treat_numeric_only::Bool)
#         gates_t = Dict(qubits_t[1] => (I_Gate{T}, X_Gate{T}))
#         gates_c = Dict(qubits_c[1] => (_00_Gate{T}, _11_Gate{T}))
#         println(nameof(I_Gate))
        
#         base_gate = make_MultiQubitQuantumGate(name=name, name_short="CX",qubits_t=qubits_t, qubits_c=qubits_c, step=step, is_treat_numeric_only=is_treat_numeric_only,
#         parameters=nothing, gates_t=gates_t, gates_c=gates_c)
#         return base_gate
#     end
# end