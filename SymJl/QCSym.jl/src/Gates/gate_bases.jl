import Symbolics
import SymbolicUtils
import ..BitsRegs.AbstractBit
import ..BitsRegs.MapBitID
import ..BitsRegs.Bit

abstract type AbstractGate end
abstract type AbstractQuantumGate{T} <: AbstractGate end
abstract type AbstractSingleQubitQuantumGate{T} <: AbstractQuantumGate{T} end
abstract type AbstractMultiQubitQuantumGate{T} <: AbstractQuantumGate{T} end
abstract type AbstractInternalSingleQubitQuantumGate{T} <: AbstractSingleQubitQuantumGate{T} end

macro insert_fields_AbstractQuantumGate()
    quote
        $(esc(:(num_qubits::Int)))
        $(esc(:(num_qubits_t::Int)))
        $(esc(:(num_qubits_c::Union{Nothing, Int})))
        $(esc(:(is_parametric::Bool)))
        $(esc(:(is_treat_numeric_only::Bool)))
        $(esc(:(is_treat_alt_only::Bool)))
        $(esc(:(name::String)))
        $(esc(:(name_short::String)))
        $(esc(:(shape::Tuple{Int, Int})))
        $(esc(:(qubits::Array{T,1})))
        $(esc(:(qubits_t::Array{T,1})))
        $(esc(:(qubits_c::Union{Nothing, Array{T,1}})))
        $(esc(:(step::Int)))
        $(esc(:(num_summands_decomposed::Int)))
        $(esc(:(parameters::Union{Nothing, Dict{Symbolics.Num, <:Real}})))
        $(esc(:(atomics::Vector{<:Symbolics.Num})))
        $(esc(:(atomics_alt::Union{Nothing, Vector{<:Symbolics.Num}})))
        $(esc(:(matrix::Symbolics.Arr{Symbolics.Num,2})))
        $(esc(:(matrix_alt::Union{Nothing, Matrix{Symbolics.Num}, Matrix{Complex{Symbolics.Num}}, Symbolics.Arr{Symbolics.Num,2}, SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymbolicUtils.SymReal}})))
        $(esc(:(ids_matrix_zeros::Union{Nothing, Array{Int, 2}})))
        $(esc(:(matrix_numeric::Union{Nothing, Array{Complex,2}})))
        $(esc(:(matrix22_t::Union{Nothing, Dict{Int, Vector{Symbolics.Arr{Symbolics.Num,2}}}})))
        $(esc(:(matrix22_t_alt::Union{Nothing, Dict{Int, Vector{Symbolics.Arr{Symbolics.Num,2}}}})))
        $(esc(:(matrix22_c::Union{Nothing, Dict{Int, Vector{Symbolics.Arr{Symbolics.Num,2}}}})))
        $(esc(:(matrix22_t_numeric::Union{Nothing, Dict{Int, Vector{Array{Complex,2}}}})))
        $(esc(:(matrix22_c_numeric::Union{Nothing, Dict{Int, Vector{Array{Complex,2}}}})))
    end
end

mutable struct mutable_BaseQuantumGate_for_construction{T<:AbstractBit} <: AbstractQuantumGate{T}
    @insert_fields_AbstractQuantumGate()
    
    # function mutable_BaseQuantumGate_for_construction(;name="", name_short="", shape=(0,0),
    # qubits_t=[Bit(context=MapBitID(), index_local=-1, is_quantum=true, name_reg="default")], qubits_c=nothing,
    # step=0, num_summands_decomposed=0, parameters=nothing, is_treat_numeric_only=false,
    # ids_matrix_zeros=nothing, matrix_numeric=nothing, matrix22_t=nothing, matrix22_t_alt=nothing,
    # matrix22_c=nothing, matrix22_t_numeric=nothing, matrix22_c_numeric=nothing)
    function mutable_BaseQuantumGate_for_construction(; 
        is_treat_numeric_only::Bool,
        name_prefix::String, name_short::String, qubits_t::Vector{T}, qubits_c::Union{Nothing, Vector{T}}=nothing,
        step::Int, num_summands_decomposed::Int, parameters::Union{Nothing, Vector{Symbolics.Num}, Dict{Symbolics.Num, <:Complex}, Dict{Symbolics.Num, <:Real}}=nothing) where {T<:AbstractBit}
        
        num_qubits, num_qubits_t, num_qubits_c = _determine_num_qubits(qubits_t, qubits_c)
        is_parametric = parameters !== nothing
        is_treat_numeric_only = is_treat_numeric_only
        is_treat_alt_only = false
        name=_generate_name_str(name_prefix*name_short, step, qubits_t, qubits_c)
        name_short=name_short
        shape=(2^num_qubits, 2^num_qubits)
        qubits_t=qubits_t
        qubits_c=qubits_c
        qubits = isnothing(qubits_c) ? qubits_t : vcat(qubits_t, qubits_c)
        step=step
        num_summands_decomposed=num_summands_decomposed
        if is_parametric
            if parameters isa Vector{Symbolics.Num}
                parameters = Dict(p => 0.0 for p in parameters)
            elseif parameters isa Dict{Symbolics.Num, <:Complex}
                parameters = parameters
            elseif parameters isa Dict{Symbolics.Num, <:Real}
                parameters = parameters
            else
                error("parameters must be either a Vector{Symbolics.Num} or a Dict{Symbolics.Num, <:Complex} or a Dict{Symbolics.Num, <:Real}")
            end
        end        
        atomics = nothing # will be set after matrix
        atomics_alt = parameters === nothing ? nothing : collect(keys(parameters))
        matrix = eval(Meta.parse("Symbolics.@variables($(name)[1:$(shape[1]),1:$(shape[2])])"))
        matrix = matrix[1]
        matrix_alt = nothing
        atomics = Symbolics.scalarize(matrix)[:]
        
        ids_matrix_zeros=nothing
        matrix_numeric=nothing
        matrix22_t=nothing
        matrix22_t_alt=nothing
        matrix22_c=nothing
        matrix22_t_numeric=nothing
        matrix22_c_numeric=nothing

        

        return new{T}(num_qubits, num_qubits_t, num_qubits_c,
            is_parametric, is_treat_numeric_only, is_treat_alt_only,
            name, name_short, shape, qubits, qubits_t, qubits_c,
            step, num_summands_decomposed, parameters, atomics, atomics_alt,
            matrix, matrix_alt, ids_matrix_zeros, matrix_numeric,
            matrix22_t, matrix22_t_alt, matrix22_c,
            matrix22_t_numeric, matrix22_c_numeric
        )
    end
end


"""
    @constructor_from_mutable_base(struct_name, base_varname)

Generates an inner constructor for an immutable struct named `struct_name`.
The generated constructor accepts a single `mutable_BaseQuantumGate_for_construction`
instance (bound to parameter `base_varname`) and copies all its fields positionally
into `new(...)`.

Intended to be placed inside a struct body, where a separate user-facing constructor
builds and mutates a `mutable_BaseQuantumGate_for_construction`, then delegates to
this generated constructor.

# Usage

    struct MyGate{T<:AbstractBit} <: AbstractQuantumGate{T}
        @insert_fields_AbstractQuantumGate()
        @constructor_from_mutable_base(MyGate, base_gate)

        function MyGate(; qubits_t, qubits_c=nothing, step)
            base_gate = mutable_BaseQuantumGate_for_construction(...)
            # ... mutate base_gate fields as needed ...
            return MyGate(base_gate)   # dispatches to the generated constructor
        end
    end
"""
macro constructor_from_mutable_base(struct_name, base_varname)
    quote
        function $(esc(struct_name))($(esc(base_varname))::mutable_BaseQuantumGate_for_construction{_COMB_T}) where {_COMB_T <: AbstractBit}
            $(esc(:new)){_COMB_T}(Tuple(getfield($(esc(base_varname)), f) for f in fieldnames(typeof($(esc(base_varname)))))...)
        end
    end
end


function _generate_name_str(name::AbstractString, step::Int, qubits_t::Array{T,1}, qubits_c::Union{Nothing, Array{T,1}}) where {T<:AbstractBit}
    return name * "_s" * string(step) * "qt" * join((q.index_global for q in qubits_t), "") * "qc" * (isnothing(qubits_c) ? "" : join((q.index_global for q in qubits_c), ""))
end

function _determine_num_qubits(qubits_t::Array{T,1}, qubits_c::Union{Nothing, Array{T,1}}) where {T<:AbstractBit}
    num_qubits_t = size(qubits_t,1)
    num_qubits_c = isnothing(qubits_c) ? 0 : size(qubits_c,1)
    return num_qubits_t + num_qubits_c, num_qubits_t, num_qubits_c
end



#@kwdef struct BaseQuantumGate{T<:AbstractBit} <: AbstractQuantumGate{T}
#    @insert_fields_AbstractQuantumGate()
#end

#println(fieldnames(BaseQuantumGate))

# @kwdef struct BaseQuantumGate{T<:AbstractBit} <: AbstractQuantumGate
#     #is_should_be_listed_in_gate_collection::Bool
#     #name_gate_collection::Union{Nothing, String}
#     num_qubits::Int
#     num_qubits_t::Int
#     num_qubits_c::Union{Nothing, Int}
#     is_parametric::Bool
#     is_treat_numeric_only::Bool
    
#     name::String
#     name_short::String
#     shape::Tuple{Int, Int}
#     qubits::                 Array{T,1}
#     qubits_t::               Array{T,1}
#     qubits_c::Union{Nothing, Array{T,1}}
#     step::Int
#     num_summands_decomposed::Int
#     parameters::Union{Nothing, Vector{Dict{String, Complex}}}
#     atomics::                   Vector{<:Symbolics.Num}
#     atomics_alt::Union{Nothing, Vector{<:Symbolics.Num}}
#     matrix::                   Symbolics.Arr{Symbolics.Num,2}
#     matrix_alt::Union{Nothing, Matrix{Symbolics.Num}, Symbolics.Arr{Symbolics.Num,2}, SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymbolicUtils.SymReal}}
#     ids_matrix_zeros::Union{Nothing, Array{Int, 2}}
#     matrix_numeric::Union{Nothing, Array{Complex,2}}
#     matrix22_t::    Union{Nothing, Dict{Int, Vector{Symbolics.Arr{Symbolics.Num,2}}}}
#     matrix22_t_alt::Union{Nothing, Dict{Int, Vector{Symbolics.Arr{Symbolics.Num,2}}}}
#     matrix22_c::    Union{Nothing, Dict{Int, Vector{Symbolics.Arr{Symbolics.Num,2}}}}
#     matrix22_t_numeric::Union{Nothing, Dict{Int, Vector{Array{Complex,2}}}}
#     matrix22_c_numeric::Union{Nothing, Dict{Int, Vector{Array{Complex,2}}}}
# end


Base.show(io::IO, gate::T) where {T <: AbstractQuantumGate} = begin
    println(io, "Quantum Gate: ", gate.name)
    for name in fieldnames(typeof(gate))
        println(io, "  ", name, ": ", getfield(gate, name))
    end
end
