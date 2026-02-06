import Symbolics
import SymbolicUtils

abstract type AbstractGate end
abstract type AbstractQuantumGate <: AbstractGate end
abstract type AbstractSingleQubitQuantumGate <: AbstractQuantumGate end
abstract type AbstractMultiQubitQuantumGate <: AbstractQuantumGate end
abstract type AbstractInternalSingleQubitQuantumGate <: AbstractSingleQubitQuantumGate end

@kwdef mutable struct BaseQuantumGate <: AbstractQuantumGate
    #is_should_be_listed_in_gate_collection::Bool
    #name_gate_collection::Union{Nothing, String}
    num_qubits::UInt
    num_qubits_t::UInt
    num_qubits_c::Union{Nothing, UInt}
    is_parametric::Bool
    is_treat_numeric_only::Bool
    
    name::String
    name_short::String
    shape::Tuple{Int, Int}
    qubits::                 Array{Int,1}
    qubits_t::               Array{Int,1}
    qubits_c::Union{Nothing, Array{Int,1}}
    step::UInt
    num_summands_decomposed::UInt
    parameters::Union{Nothing, Vector{Dict{String, Complex}}}
    atomics::                   Vector{<:Symbolics.Num}
    atomics_alt::Union{Nothing, Vector{<:Symbolics.Num}}
    matrix::                   Symbolics.Arr{Symbolics.Num,2}
    matrix_alt::Union{Nothing, Matrix{Symbolics.Num}, Symbolics.Arr{Symbolics.Num,2}, SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymbolicUtils.SymReal}}
    ids_matrix_zeros::Union{Nothing, Array{Int, 2}}
    matrix_numeric::Union{Nothing, Array{Complex,2}}
    matrix22_t::    Union{Nothing, Dict{UInt, Vector{Symbolics.Arr{Symbolics.Num,2}}}}
    matrix22_t_alt::Union{Nothing, Dict{UInt, Vector{Symbolics.Arr{Symbolics.Num,2}}}}
    matrix22_c::    Union{Nothing, Dict{UInt, Vector{Symbolics.Arr{Symbolics.Num,2}}}}
    matrix22_t_numeric::Union{Nothing, Dict{UInt, Vector{Array{Complex,2}}}}
    matrix22_c_numeric::Union{Nothing, Dict{UInt, Vector{Array{Complex,2}}}}
end

function _generate_name_str(name::AbstractString, step::UInt, qubits_t::Array{Int,1}, qubits_c::Union{Nothing, Array{Int,1}})
    return name * "_s" * string(step) * "qt" * join(qubits_t, "") * "qc" * (isnothing(qubits_c) ? "" : join(qubits_c, ""))
end

function _determine_num_qubits(qubits_t::Array{Int,1}, qubits_c::Union{Nothing, Array{Int,1}})
    num_qubits_t = size(qubits_t,1)
    num_qubits_c = isnothing(qubits_c) ? 0 : size(qubits_c,1)
    return num_qubits_t + num_qubits_c, num_qubits_t, num_qubits_c
end

function make_BaseQuantumGate(;name="", name_short="", shape=(0,0),
    qubits_t=zeros(Int,0), qubits_c=nothing,
    step=0, num_summands_decomposed=0, parameters=nothing, is_treat_numeric_only=false,
    ids_matrix_zeros=nothing, matrix_numeric=nothing, matrix22_t=nothing, matrix22_t_alt=nothing,
    matrix22_c=nothing, matrix22_t_numeric=nothing, matrix22_c_numeric=nothing
    )
    tmp_name = _generate_name_str(name, step, qubits_t, qubits_c)
    num_qubits, num_qubits_t, num_qubits_c = _determine_num_qubits(qubits_t, qubits_c)
    is_parametric = !isnothing(parameters)
    qubits = vcat(qubits_t, isnothing(qubits_c) ? Int[] : qubits_c)
    atomics = nothing # will be set after matrix
    atomics_alt = nothing
    matrix = eval(Meta.parse("Symbolics.@variables($(tmp_name)[1:$(shape[1]),1:$(shape[2])])"))
    matrix = matrix[1]
    matrix_alt = nothing
    atomics = Symbolics.scalarize(matrix)[:]
    
    
    return BaseQuantumGate(num_qubits=num_qubits, num_qubits_t=num_qubits_t, num_qubits_c=num_qubits_c,
        is_parametric=is_parametric, is_treat_numeric_only=is_treat_numeric_only,
        name=name, name_short=name_short, shape=shape, qubits=qubits, qubits_t=qubits_t, qubits_c=qubits_c,
        step=step, num_summands_decomposed=num_summands_decomposed, parameters=parameters, atomics=atomics, atomics_alt=atomics_alt,
        matrix=matrix, matrix_alt=matrix_alt, ids_matrix_zeros=ids_matrix_zeros, matrix_numeric=matrix_numeric,
        matrix22_t=matrix22_t, matrix22_t_alt=matrix22_t_alt, matrix22_c=matrix22_c,
        matrix22_t_numeric=matrix22_t_numeric, matrix22_c_numeric=matrix22_c_numeric
    )
end

function make_SingleQubitQuantumGate(;name="", name_short="",
    qubits_t=0, qubits_c=nothing,
    step=0, num_summands_decomposed=0, parameters=nothing, is_treat_numeric_only=false,
    ids_matrix_zeros=nothing, matrix_numeric=nothing, matrix22_t=nothing, matrix22_t_alt=nothing,
    matrix22_c=nothing, matrix22_t_numeric=nothing, matrix22_c_numeric=nothing
    )
    qubits_t = Array{Uint,1}(qubits_t)
    base_gate = make_BaseQuantumGate(name=name, name_short=name_short, shape=(2,2),
    qubits_t=qubits_t, qubits_c=qubits_c,
    step=step, num_summands_decomposed=num_summands_decomposed, parameters=parameters, is_treat_numeric_only=is_treat_numeric_only,
    ids_matrix_zeros=ids_matrix_zeros, matrix_numeric=matrix_numeric, matrix22_t=matrix22_t, matrix22_t_alt=matrix22_t_alt,
    matrix22_c=matrix22_c, matrix22_t_numeric=matrix22_t_numeric, matrix22_c_numeric=matrix22_c_numeric
    )
    base_gate.matrix22_t = Dict{qubits_t[1] => base_gate.matrix}
    base_gate.matrix22_t_numeric = Dict{qubits_t[1] => base_gate.matrix_numeric}
    return base_gate
end

function make_MultiQubitQuantumGate(;name::String, name_short::String, qubits_t::Array{Int,1}, qubits_c::Array{Int,1}, step::UInt, is_treat_numeric_only::Bool, parameters::Union{Nothing, Vector{Dict{String, Complex}}}, gates_t, gates_c, ids_matrix_zeros=nothing, matrix_numeric=nothing)
    tmp_name = _generate_name_str(name, step, qubits_t, qubits_c)
    println(tmp_name)
    num_qubits, num_qubits_t, num_qubits_c = _determine_num_qubits(qubits_t, qubits_c)
    shape = (2^num_qubits, 2^num_qubits)
    is_parametric = !isnothing(parameters)
    num_summands_decomposed = length(gates_t[qubits_t[1]])
    qubits = vcat(qubits_t, qubits_c)
    atomics = nothing # will be set after matrix
    atomics_alt = nothing
    matrix = eval(Meta.parse("Symbolics.@variables($(tmp_name)[1:$(2^num_qubits),1:$(2^num_qubits)])"))
    matrix = matrix[1]
    matrix_alt = nothing
    atomics = Symbolics.scalarize(matrix)[:]
    # matrix_numeric
    #tmp_gates22_t = Dict(qt => collect(gtt(string(nameof(gtt))*"_"*tmp_name, Array{UInt,1}(gt),step,is_treat_numeric_only) for gtt in gt) for (qt, gt) in gates_t)
    #tmp_gates22_c = Dict(qc => collect(gcc(string(nameof(gcc))*"_"*tmp_name, Array{UInt,1}([gc]),step,is_treat_numeric_only) for gcc in gc) for (qc, gc) in gates_c)tmp_gates22_t = Dict()
    tmp_gates22_t = Dict()
    tmp_gates22_c = Dict()
    for (qt, gt) in gates_t
        tmp_gates22_t[qt] = []
        for gtt in gt
            _name = first(rsplit(replace(string(nameof(gtt)), "_Gate" => "") * "_" * tmp_name, "_", limit=2))
            push!(tmp_gates22_t[qt], gtt(name=_name, qubits_t=Array{Int,1}([qt]), step=step, is_treat_numeric_only=is_treat_numeric_only))
        end
    end
    for (qc, gc) in gates_c
        tmp_gates22_c[qc] = []
        for gcc in gc
            _name = first(rsplit(replace(string(nameof(gcc)), "_Gate" => "") * "_" * tmp_name, "_", limit=2))
            push!(tmp_gates22_c[qc], gcc(name=_name, qubits_t=Array{Int,1}([qc]), step=step, is_treat_numeric_only=is_treat_numeric_only))
        end
    end
    matrix22_t = Dict(q => collect(getfield(gg, :matrix) for gg in g) for (q, g) in tmp_gates22_t)
    matrix22_c = Dict(q => collect(getfield(gg, :matrix) for gg in g) for (q, g) in tmp_gates22_c)
    matrix22_t_numeric = Dict(q => collect(getfield(gg, :matrix_numeric) for gg in g) for (q, g) in tmp_gates22_t)
    matrix22_c_numeric = Dict(q => collect(getfield(gg, :matrix_numeric) for gg in g) for (q, g) in tmp_gates22_c)
    
    matrix22_t_alt = nothing
    return BaseQuantumGate(num_qubits=num_qubits, num_qubits_t=num_qubits_t, num_qubits_c=num_qubits_c,
        is_parametric=is_parametric, is_treat_numeric_only=is_treat_numeric_only,
        name=name, name_short=name_short, shape=shape, qubits=qubits, qubits_t=qubits_t, qubits_c=qubits_c,
        step=step, num_summands_decomposed=num_summands_decomposed, parameters=parameters, atomics=atomics, atomics_alt=atomics_alt,
        matrix=matrix, matrix_alt=matrix_alt, ids_matrix_zeros=ids_matrix_zeros, matrix_numeric=matrix_numeric,
        matrix22_t=matrix22_t, matrix22_t_alt=matrix22_t_alt, matrix22_c=matrix22_c,
        matrix22_t_numeric=matrix22_t_numeric, matrix22_c_numeric=matrix22_c_numeric
    )
end

Base.show(io::IO, gate::BaseQuantumGate) = begin
    println(io, "Quantum Gate: ", gate.name)
    for name in fieldnames(typeof(gate))
        println(io, "  ", name, ": ", getfield(gate, name))
    end    
end
