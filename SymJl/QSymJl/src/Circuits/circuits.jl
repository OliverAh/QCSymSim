import Symbolics
import ..BitsRegs
import ..Gates

struct GateCollection
    collections::Dict{typeof(Gates.AbstractGate), Vector{Any}}
end

@kwdef struct QuantumCircuit
    name::String = "MyCircuit"
    context::BitsRegs.MapBitID = BitsRegs.MapBitID()
    qregs::Vector{BitsRegs.QReg} = BitsRegs.QReg[]
    cregs::Vector{BitsRegs.CReg} = BitsRegs.CReg[]
    gatecollection::GateCollection = GateCollection(Dict())
    barriers::Vector{Any} = Any[]
end

function get_num_qubits(qc::QuantumCircuit)
    return sum(q.size for q in qc.qregs)
end

function add_qreg(circuit::QuantumCircuit, name::String="", size::Int=99)
    @assert !has_qreg(circuit.context["q"]["l2g"], name) "Quantum register with name $name already exists"
    qreg = BitsRegs.QReg(context=circuit.context, name=name, size=size)
    push!(circuit.qregs, qreg)
    return qreg
end

function add_creg(circuit::QuantumCircuit, name::String="", size::Int=99)
    @assert !has_creg(circuit.context["c"]["l2g"], name) "Classical register with name $name already exists"
    creg = BitsRegs.CReg(context=circuit.context, name=name, size=size)
    push!(circuit.cregs, creg)
    return creg
end

"""kwargs are exclusively forwarded to the added gate. See their respective documentation."""
function add_gate(circuit::QuantumCircuit, gate::typeof(Gates.AbstractGate); kwargs...)
    if !haskey(circuit.gatecollection.collections, gate)
        circuit.gatecollection.collections[gate] = Any[]
    end
    push!(circuit.gatecollection.collections[gate], gate(; kwargs...))
end

function assemble_symbolic_unitary(qc::QuantumCircuit, replace_symbolic_zeros::Bool, replace_symbolic_ones::Bool)
    U = Symbolics.@variables U[1:2^get_num_qubits(qc), 1:2^get_num_qubits(qc)]

end